//
//
// FLIRTLib - Fast Laser Interesting Region Transform Library
// Copyright (C) 2009-2010 Gian Diego Tipaldi and Kai O. Arras
//
// This file is part of FLIRTLib.
//
// FLIRTLib is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// FLIRTLib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with FLIRTLib.  If not, see <http://www.gnu.org/licenses/>.
//

#include "CurvatureDetector.h"

#include <utils/Regression.h>


#include <boost/version.hpp>

#if ((BOOST_VERSION / 100) % 1000) == 54
// Workaround to use dijkstra_shortest_paths which has a bug in Boost 1.54 from Boost 1.55
// dijkstra_shortest_paths is used prim_minimum_spanning_tree.hpp
// Here we use the 1.55 version instead (downloaded into the current folder):
//
// wget https://raw.githubusercontent.com/boostorg/graph/boost-1.55.0/include/boost/graph/dijkstra_shortest_paths.hpp
// Including the local dijkstra_shortest_paths here will shadow any subsequent includes of the system 1.54 version
#warning "Using dijkstra_shortest_paths.hpp from Boost 1.55!"
#include "dijkstra_shortest_paths.hpp"
#endif

#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>

CurvatureDetector::CurvatureDetector(const PeakFinder* peak, unsigned int scales, double sigma, double step, unsigned int dmst):
    m_peakFinder(peak),
    m_scaleNumber(scales),
    m_baseSigma(sigma),
    m_sigmaStep(step),
    m_useMaxRange(false),
    m_dmst(dmst)
{
    //计算尺度库m_scales
    computeScaleBank();
}

//detect重载1
unsigned int CurvatureDetector::detect(const LaserReading& reading, std::vector<InterestPoint*>& _point) const
{
    Graph graph;

    //operatorA为 scale数量×顶点数量 的二维数组
    //里面存放各个尺度下平滑后的结果
    std::vector< std::vector<Point2D> > operatorA;

    //signalDiff为 scale数量×顶点数量 的二维数组
    //里面存放平滑后的点到该点原位置的距离
    std::vector< std::vector<double> > signalDiff;

    //indexes存放候选的关键点
    std::vector< std::vector<unsigned int> > indexes;

    //调用detect重载2
    return detect(reading, _point, graph, operatorA, signalDiff, indexes);
}

//detect重载2
unsigned int CurvatureDetector::detect(const LaserReading& reading, std::vector<InterestPoint*>& point,
				       Graph& graph,
				       std::vector< std::vector<Point2D> >& operatorA, 
				       std::vector< std::vector<double> >& signalDiff,
				       std::vector< std::vector<unsigned int> >& indexes) const
{
    //maxRangeMapping提供从图顶点下标到激光点下标的映射，图中第k个顶点在激光中为第maxRangeMapping[k]个点
    //图顶点数量可能少于激光点数量，因为图顶点可能不包含最大射程的激光点
    std::vector<unsigned int> maxRangeMapping;

    //graphPoints用来存放所有激光点的世界坐标
    std::vector<Point2D> graphPoints;

    //计算图
    //reading为输入参数
    //graphPoints、graph、maxRangeMapping为输出参数
    computeGraph(reading, graphPoints, graph, maxRangeMapping);

    //调用detect重载4
    //graph、graphPoints作为输入参数
    //operatorA、signalDiff、indexs作为输出参数
    detect(graph, graphPoints, operatorA, signalDiff, indexes);

    return computeInterestPoints(reading, operatorA, point, indexes, maxRangeMapping);
}

//detect重载3
unsigned int CurvatureDetector::detect( const LaserReading& reading, std::vector<InterestPoint*>& point,
					std::vector< double >& signal,
					std::vector< std::vector<double> >& signalSmooth,
					std::vector< std::vector<double> >& signalDiff,
					std::vector< std::vector<unsigned int> >& indexes) const
{
    std::vector<unsigned int> maxRangeMapping;
    std::vector<Point2D> graphPoints;
    Graph graph;
    std::vector< std::vector<Point2D> > operatorA;
    computeGraph(reading, graphPoints, graph, maxRangeMapping);
    //调用detect重载4
    detect(graph, graphPoints, operatorA, signalDiff, indexes);
    signal.resize(graphPoints.size());
    signalSmooth.resize(m_scales.size(), std::vector<double> (graphPoints.size()));
    for(unsigned int i = 0; i < graphPoints.size(); i++){
		Point2D difference = graphPoints[i] - reading.getLaserPose();
		signal[i] = hypot(difference.x, difference.y);
		for(unsigned int scale = 0; scale < m_scales.size(); scale++){
			difference = operatorA[scale][i] - reading.getLaserPose();
			signalSmooth[scale][i] = hypot(difference.x, difference.y);
		}
    }
    return computeInterestPoints(reading, operatorA, point, indexes, maxRangeMapping);
}

void CurvatureDetector::computeGraph(const LaserReading& reading, std::vector<Point2D>& graphPoints, Graph& graph, std::vector<unsigned int>& maxRangeMapping) const
{
    //获取所有激光点的世界坐标
    const std::vector<Point2D>& worldPoints = reading.getWorldCartesian();

    //graphPoints和worldPoints相比，可能少了最大量程的那些点
    graphPoints.reserve(worldPoints.size());

    //创建临时存放图的“边”及其“权值”的变量
    std::vector<GraphEdge> edges;
    std::vector< boost::property < boost::edge_weight_t, double > > weights;
    edges.reserve(worldPoints.size()*worldPoints.size()); //预留充足的容量（capacity）
    weights.reserve(worldPoints.size()*worldPoints.size());

    unsigned int currentVertexNumber = 0; //当前的顶点编号
    //对所有激光点进行遍历
    //遍历完的结果：edges中存放了所有顶点之间的连接边，weights中存放相应的边的权值
    for(unsigned int i = 0; i < worldPoints.size(); i++){
        //如果m_useMaxRange=false，则忽略最大量程的激光点
        if(m_useMaxRange || reading.getRho()[i] < reading.getMaxRange()){
            //将当前激光点加入到图的顶点数组
            graphPoints.push_back(worldPoints[i]);

            maxRangeMapping.push_back(i); //如果m_useMaxRange为false，则maxRangeMapping不包含最大测量的下标

            unsigned int targetVertexNumber = currentVertexNumber + 1;

            //建立当前遍历的第i个点与它后面所有点之间的连接边
            for(unsigned int j = i + 1; j < worldPoints.size(); j++){
                if(m_useMaxRange || reading.getRho()[j] < reading.getMaxRange()){
                    //边的权值即为两个激光点之间的欧氏距离
                    Point2D difference = worldPoints[i] - worldPoints[j];
                    double weight = hypot(difference.x, difference.y);

                    //将当前生成的边和其权值分别存入相应的数组
                    //GraphEdge(a, b)是使用std::pair()的构造函数，a、b表示这条边连接的两个顶点
                    edges.push_back(GraphEdge(currentVertexNumber,targetVertexNumber));
                    weights.push_back(weight);
                    targetVertexNumber++;
                }
            }
            currentVertexNumber++;
        }
    }
    
    //costGraph是一个邻接矩阵形式存储的图，用它重新表示上文以graphPoints、edges、weights三个数组表示的图
    MatrixGraph costGraph(currentVertexNumber); //参数为图顶点的数量
    for(unsigned int i = 0; i < edges.size(); i++){
        //向costGraph中添加边
        boost::add_edge(edges[i].first, edges[i].second, weights[i], costGraph);
    }
    boost::remove_edge(0, currentVertexNumber - 1, costGraph); //移除第一个顶点和最后一个顶点之间的边，作用未知
    
    //以多次MST操作的方式构建DMST
    for(unsigned int iter = 0; iter <= m_dmst; iter++){
        //根据prim_minimum_spanning_tree()输出参数的需要，创建顶点数组，用来保存最小生成树的结果
        //predecessor是“先驱”的意思，保存着节点的父节点下标，(predecessor[u], u)为一个边，如果predecessor[r]=r，说明r为根节点
        std::vector < boost::graph_traits < Graph >::vertex_descriptor > predecessor(boost::num_vertices(costGraph)); //num_vertices()用于获取一个图的顶点数

        //使用Prim算法求得最小生成树
        boost::prim_minimum_spanning_tree(costGraph, &predecessor[0]); //第二个参数的这种写法是规范写法

        //对最小生成树的所有边进行遍历，默认将第0个点视作树根，对(predecessor[0], 0)不遍历
        //遍历的目的是：① 将最小生成树的边添加到图graph中；②将最生成树的边从图costGraph中移除
        for(unsigned int i = 1; i < predecessor.size(); i++){
            boost::add_edge(predecessor[i], i, boost::get(boost::edge_weight, costGraph, boost::edge(predecessor[i], i, costGraph).first), graph);
            boost::remove_edge(predecessor[i], i, costGraph);
        }
    }
}

//detect重载4
/*
输入参数
  graph：DMST图
  graphPoints：图中所有顶点的世界坐标
输出参数
  operatorA：在每个尺度下曲线平滑后的结果
  signalDiff：在每个尺度下每个点平滑后与平滑前的距离差
  indexes：在每个尺度下候选关键点的下标
*/
void CurvatureDetector::detect(const Graph& graph, const std::vector<Point2D>& graphPoints, std::vector< std::vector<Point2D> >& operatorA, std::vector< std::vector<double> >& signalDiff, 
				       std::vector< std::vector<unsigned int> >& indexes) const
{
    operatorA.resize(m_scales.size());
    signalDiff.resize(m_scales.size());
    indexes.resize(m_scales.size());
    
    //获取图的顶点数量
    unsigned int vertexNumber = boost::num_vertices(graph);

    //使用Johnson全源最短路径算法 获取所有顶点之间的通行代价 保存到distances二维数组中
    std::vector< std::vector< double > > distances(vertexNumber, std::vector<double>(vertexNumber));
    boost::johnson_all_pairs_shortest_paths(graph, distances);
   
    //对尺度库中的所有尺度进行遍历
    for(unsigned int scale = 0; scale < m_scales.size(); scale++){
        double currentScale = m_scales[scale]; //即正态分布的标准差σ
        double normalizer = sqrt(2*M_PI)*currentScale; //即高斯函数中e^(...)前面的系数

        std::vector<double> densities(vertexNumber, 0.); //用来存放每个点周围的疏密程度

        operatorA[scale].resize(vertexNumber, Point2D());
        signalDiff[scale].resize(vertexNumber,0.);

        double weights[vertexNumber][vertexNumber];
        double weightNormalizer[vertexNumber];

        for(unsigned int vertex1 = 0; vertex1 < vertexNumber; vertex1++){
            for(unsigned int vertex2 = 0; vertex2 < vertexNumber; vertex2++){
                //利用高斯函数求vertex1点周围的疏密程度
                weights[vertex1][vertex2] = normalizer * exp(-distances[vertex1][vertex2] * distances[vertex1][vertex2]/(2 * currentScale * currentScale));
                densities[vertex1] += weights[vertex1][vertex2];
            }
        }

        for(unsigned int vertex1 = 0; vertex1 < vertexNumber; vertex1++){
            weightNormalizer[vertex1] = 0.;
            for(unsigned int vertex2 = 0; vertex2 < vertexNumber; vertex2++){
                //将“vertex2对vertex1贡献的权值”除以“vertex1与vertex2密度之积”完成权值的标准化
                weights[vertex1][vertex2] /= densities[vertex1] * densities[vertex2];
                weightNormalizer[vertex1] += weights[vertex1][vertex2];
            }
        }

        for(unsigned int vertex1 = 0; vertex1 < vertexNumber; vertex1++){
            //对vertex1的(x, y)坐标进行高斯平滑
            for(unsigned int vertex2 = 0; vertex2 < vertexNumber; vertex2++){
                operatorA[scale][vertex1] = operatorA[scale][vertex1] + weights[vertex1][vertex2] * graphPoints[vertex2];
            }
            //（这一步不理解）
            operatorA[scale][vertex1] = operatorA[scale][vertex1] * ( 1. / weightNormalizer[vertex1]);

            Point2D pointDifference = operatorA[scale][vertex1] - graphPoints[vertex1];
            double temporary = 2 * (1/currentScale) * hypot(pointDifference.x, pointDifference.y);
            signalDiff[scale][vertex1] =  temporary * exp(-temporary); //（为什么使用x*e^(-x)函数？）
        }

        for(unsigned int j = 2; j < signalDiff[scale].size() - 2; j++){ //不包括scan最边缘的两个点
            if(m_peakFinder->isPeak(signalDiff[scale], j)){
                indexes[scale].push_back(j);
            }
        }
    }
}

unsigned int CurvatureDetector::computeInterestPoints(const LaserReading& reading, const std::vector< std::vector<Point2D> >& operatorA, std::vector<InterestPoint*>& point, 
                                                      const std::vector< std::vector<unsigned int> >& indexes, std::vector<unsigned int>& maxRangeMapping) const
{
    point.clear();
    point.reserve(reading.getRho().size());
    const std::vector<Point2D>& worldPoints = reading.getWorldCartesian();
    for(unsigned int i = 0; i < indexes.size(); i++){ //i表示第几个尺度(scale)
        for(unsigned int j = 0; j < indexes[i].size(); j++){
            OrientedPoint2D pose;
            unsigned int pointIndex = maxRangeMapping[indexes[i][j]]; //maxRangeMapping提供了图顶点下标到worldPoints下标的映射

            // Reomoving the detection in the background and pushing it to the foreground
            double rangeBefore = (pointIndex > 1)? reading.getRho()[pointIndex - 1] : reading.getMaxRange();
            double rangeAfter = (pointIndex < worldPoints.size() - 1)? reading.getRho()[pointIndex + 1] : reading.getMaxRange();
            double rangeCurrent = reading.getRho()[pointIndex];
            if(rangeBefore < rangeAfter){
                if(rangeBefore < rangeCurrent){
                    pointIndex = pointIndex - 1;
                }
            }
            else if(rangeAfter < rangeCurrent){
                pointIndex = pointIndex + 1;
            }

            // Removing max range reading
            if(reading.getRho()[pointIndex] >= reading.getMaxRange()){
                continue;
            }


            pose.x =  (worldPoints[pointIndex]).x;
            pose.y =  (worldPoints[pointIndex]).y;
            Point2D difference = operatorA[i][indexes[i][j]] - worldPoints[pointIndex];
            pose.theta = atan2(difference.y, difference.x);

            //非极大值抑制
            bool exists = false;
            for(unsigned int k = 0; !exists && k < point.size(); k++){
                //判断2分米内是否有其他关键点
                exists = exists || (fabs(pose.x - point[k]->getPosition().x) <= 0.2 &&  fabs(pose.y - point[k]->getPosition().y) <= 0.2);
            }
            if(exists) continue;

            double distance = 2. * m_scales[i]; //distance为2σ（2σ以内的概率为95%）
            Point2D diffStart = pose - worldPoints.front(); //std::vector的front()与bace()分别返回第一个元素和最后一个元素的引用
            Point2D diffEnd = pose - worldPoints.back();

            //如果该点与第一个点或最后一个点的距离小于2σ，则抛弃
            if(hypot(diffStart.x, diffStart.y) < distance || hypot(diffEnd.x, diffEnd.y) < distance){
                continue;
            }

            //将距离关键点2σ以内的点作为“支持点”
            std::vector<Point2D> support;
            for(unsigned int k = 0; k < worldPoints.size(); k++){
                Point2D diff2 = pose - worldPoints[k];
                if(hypot(diff2.x, diff2.y) < distance) support.push_back(worldPoints[k]);
            }

            // 	    LineParameters param = computeNormals(support);
            // 	    pose.theta = normAngle(param.alpha, - M_PI);

            InterestPoint *interest = new InterestPoint(pose, distance);
            //  InterestPoint *interest = new InterestPoint(pose, m_scales[i]);
            interest->setSupport(support); //设置关键点的支持点
            interest->setScaleLevel(i); //设置关键点的尺度等级
            point.push_back(interest);
        }
    }
    return point.size();
    
}

void CurvatureDetector::computeScaleBank()
{
    m_scales.resize(m_scaleNumber);
    for(unsigned int i = 0; i < m_scales.size(); i++){
		m_scales[i] = m_baseSigma * pow(m_sigmaStep, i);
    }
}
