#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include "geometry.h"
#include "iofile.h"

#define ITERATION_TIMES 100
#define DEBUG_PICTURE
#define DEBUG_FINDMINCUT
using namespace std;


void createDual(vector<surface>& surfaceData, vector<vector<surfacePriority> >& dual, double& avgAngd, vector<vector<surfacePriority> >& angdDual);
int dijkstra(vector<vector<double> >& dual, int s, int fnum, vector<double>& priority);
int dijkstra_prior(vector<vector<surfacePriority> >& dual, int s, int fnum, vector<double>& minDist);	//用优先级队列

void posibility(vector<vector<double>>& psb, vector<vector<double>>& dist, vector<int>& REP, int knum, int fnum, int times);
void modify(vector<vector<double>>& psb, vector<vector<double>>& dist, vector<int>& REP, int knum, int fnum, int times);

void findMinCut(vector<vector<surfacePriority>>& dualC, vector<int>& borderI, int REPJ, vector<BLOCK>& block);
void preciseDivide(vector<vector<surfacePriority>>& dualC, vector<int>& borderI, int REPI, int REPJ, vector<BLOCK>& block);
void createFlowNetwork(vector<vector<surfacePriority>>& angdDual, vector<BLOCK>& block, vector<int>& borderI, double avgAngd,
	int REPI, int REPJ);


void createDual(vector<surface>& surfaceData, vector<vector<surfacePriority> >& dual, double& avgAngd, vector<vector<surfacePriority> >& angdDual)
{
	map<side, int> sideHash;	//索引是边，键值是某个面片的编号
	int fnum = surfaceData.size();

	vector<vector<double> > dualMatrix(fnum, vector<double>(fnum, 0));	//暂时先用邻接矩阵，后面再转成邻接表

	int sideNum = 0;
	double totalGeod = 0;
	double totalAngd = 0;

	vector<vector<double> > geodMatrix(fnum, vector<double>(fnum, 0));
	vector<vector<double> > angdMatrix(fnum, vector<double>(fnum, 0));
	//先用这个矩阵储存测地距离和角距离

	for (int i = 0; i < fnum; i++)	//第一轮遍历：先记录距离总和和边的数量
	{
		//第i个面片的三条边
		side s1(surfaceData[i].vtx1, surfaceData[i].vtx2);
		side s2(surfaceData[i].vtx2, surfaceData[i].vtx3);
		side s3(surfaceData[i].vtx3, surfaceData[i].vtx1);

		pair<map<side, int>::iterator, bool> insertPair;//如果插入失败，那么这个索引处的键值就是此面的相邻面

#ifndef DEBUG_PICTURE
		cout << "surface" << i << ": " << endl;
		cout << "side1: " << s1.a << "-" << s1.b << ", side2: " << s2.a << "-" << s2.b << ", side3: " << s3.a << "-" << s3.b << endl;
#endif // !DEBUG_PICTURE

		insertPair = sideHash.insert(pair<side, int>(s1, i));
		if (!insertPair.second)
		{	//若插入失败，边s1的另一个面和此面（i）共有边s1
			int anotherSurface = sideHash[s1];

			geodMatrix[anotherSurface][i] = getGeod(surfaceData[anotherSurface], surfaceData[i]);
			geodMatrix[i][anotherSurface] = geodMatrix[anotherSurface][i];
			angdMatrix[anotherSurface][i] = getAngd(surfaceData[anotherSurface], surfaceData[i]);
			angdMatrix[i][anotherSurface] = angdMatrix[anotherSurface][i];

			totalGeod += geodMatrix[anotherSurface][i];
			totalAngd += angdMatrix[anotherSurface][i];
			sideNum++;

#ifndef DEBUG_PICTURE
			cout << "side:" << s1.a << "->" << s1.b << ", surface:" << i << "," << anotherSurface;
			cout << ", geod:" << geodMatrix[anotherSurface][i] << ", angd:" << angdMatrix[anotherSurface][i] << endl;
#endif // !DEBUG_PICTURE

		}

		insertPair = sideHash.insert(pair<side, int>(s2, i));
		if (!insertPair.second)
		{
			int anotherSurface = sideHash[s2];

			geodMatrix[anotherSurface][i] = getGeod(surfaceData[anotherSurface], surfaceData[i]);
			geodMatrix[i][anotherSurface] = geodMatrix[anotherSurface][i];
			angdMatrix[anotherSurface][i] = getAngd(surfaceData[anotherSurface], surfaceData[i]);
			angdMatrix[i][anotherSurface] = angdMatrix[anotherSurface][i];

			totalGeod += geodMatrix[anotherSurface][i];
			totalAngd += angdMatrix[anotherSurface][i];
			sideNum++;

#ifndef DEBUG_PICTURE
			cout << "side:" << s2.a << "->" << s2.b << ", surface:" << i << "," << anotherSurface;
			cout << ", geod:" << geodMatrix[anotherSurface][i] << ", angd:" << angdMatrix[anotherSurface][i] << endl;
#endif // !DEBUG_PICTURE
		}

		insertPair = sideHash.insert(pair<side, int>(s3, i));
		if (!insertPair.second)
		{
			int anotherSurface = sideHash[s3];

			geodMatrix[anotherSurface][i] = getGeod(surfaceData[anotherSurface], surfaceData[i]);
			geodMatrix[i][anotherSurface] = geodMatrix[anotherSurface][i];
			angdMatrix[anotherSurface][i] = getAngd(surfaceData[anotherSurface], surfaceData[i]);
			angdMatrix[i][anotherSurface] = angdMatrix[anotherSurface][i];

			totalGeod += geodMatrix[anotherSurface][i];
			totalAngd += angdMatrix[anotherSurface][i];
			sideNum++;

#ifndef DEBUG_PICTURE
			cout << "side:" << s3.a << "->" << s3.b << ", surface:" << i << "," << anotherSurface;
			cout << ", geod:" << geodMatrix[anotherSurface][i] << ", angd:" << angdMatrix[anotherSurface][i] << endl;
#endif // !DEBUG_PICTURE
		}
	}

	double avgGeod = totalGeod / sideNum;
	avgAngd = totalAngd / sideNum;	//计算平均角距离

	for (int i = 0; i < fnum; i++)
	{
		for (int j = 0; j < fnum; j++)
			dualMatrix[i][j] = theta * geodMatrix[i][j] / avgGeod + (1 - theta) * angdMatrix[i][j] / avgAngd;
	}

	//构建邻接表
	for (int i = 0; i < fnum; i++)
	{
		for (int j = 0; j < fnum; j++)
		{
			if (dualMatrix[i][j])
			{
				dual[i].push_back(surfacePriority(dualMatrix[i][j], j));
			}
		}
	}
	//构建表示角距离的邻接表
	for (int i = 0; i < fnum; i++)
	{
		int size = dual[i].size();
		for (int j = 0; j < size; j++)
		{
			int thisOrder = dual[i][j].order;	//指向的面片的编号
			double angd = angdMatrix[i][thisOrder];//找到对应的角距离
			angdDual[i].push_back(surfacePriority(angd, thisOrder));
		}
	}

#ifndef DEBUG_PICTURE
	cout << "sideNum: " << sideNum << endl;
	cout << "geodMatrix:" << endl;
	for (int i = 0; i < fnum; i++)
	{
		for (int j = 0; j < fnum; j++)
			cout << setw(6) << geodMatrix[i][j] << "    ";
		cout << endl;
	}
	cout << "totalGeod: " << totalGeod << ", avgGeod: " << avgGeod << endl;

	cout << "angdMatrix:" << endl;
	for (int i = 0; i < fnum; i++)
	{
		for (int j = 0; j < fnum; j++)
			cout << angdMatrix[i][j] << "    ";
		cout << endl;
	}
	cout << "totalAngd: " << totalAngd << ", avgAngd: " << avgAngd << endl;

	cout << endl;
	cout << "dualMatrix: " << endl;
	for (int i = 0; i < fnum; i++)
	{
		for (int j = 0; j < fnum; j++)
			cout << dualMatrix[i][j] << " ";
		cout << endl;
	}

	cout << endl;
	cout << "dual: " << endl;
	for (int i = 0; i < fnum; i++)
	{
		for (int j = 0; j < dual[i].size(); j++)
			cout << dual[i][j].order;
		cout << endl;
	}
#endif // !DEBUG_PICTURE


}

int dijkstra(vector<vector<double> >& dual, int s, int fnum, vector<double>& priority)
//在对偶图dual中，以s为源点，寻找最短路，fnum是矩阵的维数, priority是一个数组，储存每个点的最短路的长度
{
	vector<bool> ifVisited(fnum, false);	//用于标记已经在路径树中的节点

	priority[s] = 0;
	for (int i = 0; i < fnum; i++)
	{
		ifVisited[s] = true;	//每一轮中s是权值最小的点
		for (int j = 0; j < fnum; j++)
		{	//遍历与s相邻的点，遇到没有visit的就更新权值
			if (dual[s][j] && (!ifVisited[j]))
			{
				if (priority[j] > priority[s] + dual[s][j])
					priority[j] = priority[s] + dual[s][j];
			}
		}

		double shortest = INT_MAX;
		for (int j = 0; j < fnum; j++)
		{
			if (ifVisited[j] == false && priority[j] < shortest)
			{
				shortest = priority[j];
				s = j;
			}
		}
	}

	int farPoint = 0;
	double maxPriority = 0;
	for (int i = 0; i < fnum; i++)
	{
		if (priority[i] > maxPriority)
		{
			maxPriority = priority[i];
			farPoint = i;
		}
	}

	return farPoint;
}


int dijkstra_prior(vector<vector<surfacePriority> >& dual, int s, int fnum, vector<double>& minDist)
{//dual为原图的邻接表，priority储存源点s到所有点的最短距离,返回从s出发，最短距离最远的那个面片的编号

	priority_queue<surfacePriority> neighbor;
	minDist[s] = 0;
	vector<bool> ifVisited(fnum, false);

	neighbor.push(surfacePriority(0, s));

	while (!neighbor.empty())
	{
		surfacePriority thisSurface = neighbor.top();
		neighbor.pop();
		if (ifVisited[thisSurface.order])
			continue;	//新出队的队头可能已经访问过了，因为是可以重复入队的
		ifVisited[thisSurface.order] = true;

		int thisOrder = thisSurface.order;
		double thisPriority = thisSurface.priority;
		//当前出队的节点（即当前确认距离的节点）的信息

		minDist[thisOrder] = thisPriority;	//把每个节点的最短距离存到minDist中

		//更新这个点的邻域并放入队列
		int neighborSize = dual[thisOrder].size();
		for (int i = 0; i < neighborSize; i++)
		{
			if (!ifVisited[dual[thisOrder][i].order])
			{
				double newPriority = thisPriority + dual[thisOrder][i].priority;
				int newOrder = dual[thisOrder][i].order;
				surfacePriority newSurface(newPriority, newOrder);
				neighbor.push(newSurface);

			}
		}
	}

	double maxMinDist = 0;
	int maxOrder = 0;
	for (int i = 0; i < fnum; i++)
	{
		if (minDist[i] > maxMinDist)
		{
			maxMinDist = minDist[i];
			maxOrder = i;
		}
	}

	return maxOrder;
}



void posibility(vector<vector<double>>& psb, vector<vector<double>>& dist, vector<int>& REP, int knum, int fnum, int times)
{//传入：存储概率的矩阵，存储距离的矩阵，存储聚类中心的向量，聚类中心个数，面片个数，迭代的次数
	if (times == ITERATION_TIMES)
		return;
	for (int i = 0; i < fnum; i++)
	{//对于第i个点，计算它属于第j个中心的概率
		double sumDist = 0;
		for (int j = 0; j < knum; j++)
			sumDist += (1 / dist[i][REP[j]]); //先得到它到所有中心的总距离

		for (int j = 0; j < knum; j++)
			psb[j][i] = (1 / dist[REP[j]][i]) / sumDist;
		//第j个中心点的概率向量
	}

	for (int i = 0; i < knum; i++)
		psb[i][REP[i]] = 1;

	modify(psb, dist, REP, knum, fnum, times);

}



void modify(vector<vector<double>>& psb, vector<vector<double>>& dist, vector<int>& REP, int knum, int fnum, int times)
{
	cout << "modify: turn" << times << endl;

	vector<int> newREP(knum);
	for (int i = 0; i < knum; i++) //一次更新所有的聚类中心
	{//对于第i个聚类中心
		int thisNewREP = REP[i];
		double minSum = (double)INT_MAX;
		for (int j = 0; j < fnum; j++)
		{//对于第j个面片，每个面片进行计算到所有的面片的PD（距离*概率）之和
			double sumPD = 0.0;
			for (int k = 0; k < fnum; k++)
				sumPD += psb[i][k] * dist[k][j];

			if (sumPD < minSum)
			{	//找到聚类中心i对应的PD和最小的那个点（新的聚类中心）
				minSum = sumPD;
				thisNewREP = j;
			}
		}

		newREP[i] = thisNewREP;
	}

	vector<int> sREP = REP; //为了在排序的时候不破坏原有的REP的顺序，新复制一个一样的vector
	vector<int> sNewREP = newREP;
	sort(sREP.begin(), sREP.end());
	sort(sNewREP.begin(), sNewREP.end());

	cout << "newREP: ";
	for (int i = 0; i < knum; i++)
		cout << newREP[i] << " ";
	cout << endl;

	bool ifSame = true;
	for (int i = 0; i < knum; i++)
	{
		if (sREP[i] != sNewREP[i])
		{
			ifSame = false;
			break;
		}
	}

	if (ifSame)
		return;
	else
	{
		for (int i = 0; i < knum; i++)
			REP[i] = newREP[i];
		posibility(psb, dist, REP, knum, fnum, times + 1);
	}

}



void findMinCut(vector<vector<surfacePriority>>& dualC, vector<int>& borderI, int REPJ, vector<BLOCK>& block)
{
#ifndef DEBUG_FINDMINCUT
	cout << "a new turn" << endl;
#endif // DEBUG_FINDMINCUT



	queue<int> Q;
	int sizeI = borderI.size();
	int fnum = block.size();

	vector<int> pre(fnum, -1); //储存每个面片在路径上的前一个节点
	vector<bool> ifIn(fnum, false);

	for (int i = 0; i < sizeI; i++)
	{	//边界先入队
		Q.push(borderI[i]);
		ifIn[borderI[i]] = true;
	}

	bool ifFound = false;
	int end;
	while (!Q.empty() && !ifFound)
	{	//只要找到通路或者队列已空就结束循环
		int thisSurface = Q.front(); Q.pop();
		int neighborSize = dualC[thisSurface].size();
		for (int i = 0; i < neighborSize; i++)
		{	//看已出队的面片的邻域 
			if (dualC[thisSurface][i].priority > 0)
			{	//如果有通路
				int newSurface = dualC[thisSurface][i].order;
				if (ifIn[newSurface] == false)
				{	//如果没入过队
					pre[newSurface] = thisSurface;
					ifIn[newSurface] = true;
					if (block[newSurface].ifConfirm == true && block[newSurface].area1 == REPJ)
					{	//如果这个新的面片是终点
						ifFound = true;
						end = newSurface;
						break;
					}
					Q.push(newSurface);
				}
			}
		}
	}

	if (ifFound)
	{	//如果找到
		double minCap = INT_MAX; //路径最小水流量
		int currentSurface = end;
		while (1)
		{	//找这条路径上的最小水流量
			int preSurface = pre[currentSurface];

#ifndef DEBUG_FINDMINCUT
			cout << "->" << preSurface;
#endif // DEBUG_FINDMINCUT

			if (preSurface < 0)
				break;

			int size = dualC[preSurface].size();
			for (int i = 0; i < size; i++)
			{
				if (dualC[preSurface][i].order == currentSurface)
				{//找当前面片与前一个面片之间的边,然后看它的边权值是否是整条路上最小的
					if (dualC[preSurface][i].priority < minCap)
						minCap = dualC[preSurface][i].priority;
					break;
				}
			}
			currentSurface = preSurface;
		}

#ifndef DEBUG_FINDMINCUT
		cout << endl;
		cout << "the min cap is: " << minCap << endl;
#endif // DEBUG_FINDMINCUT

		currentSurface = end;
		while (1)
		{	//找到这条通路上的最小边权值后，所有边权值减去它
			int preSurface = pre[currentSurface];
			if (preSurface < 0)
				break;

			int size = dualC[preSurface].size();
			for (int i = 0; i < size; i++)
			{
				if (dualC[preSurface][i].order == currentSurface)//找前一个面片指向当前面片的边
				{
					dualC[preSurface][i].priority -= minCap;
					break;
				}
			}

			size = dualC[currentSurface].size();
			for (int i = 0; i < size; i++)
			{
				if (dualC[currentSurface][i].order == preSurface)//找当前面片指向前一个面片的边
				{
					dualC[currentSurface][i].priority -= minCap;
					break;
				}
			}

			currentSurface = preSurface;
		}

		findMinCut(dualC, borderI, REPJ, block);
	}
	else
	{
		return;
	}
}

void preciseDivide(vector<vector<surfacePriority>>& dualC, vector<int>& borderI, int REPI, int REPJ, vector<BLOCK>& block)
{
	queue<int> Q;

	int fnum = block.size();
	vector<bool> ifIn(fnum, false);

	int sizeI = borderI.size();
	for (int i = 0; i < sizeI; i++)
	{	//边界先入队
		Q.push(borderI[i]);
		ifIn[borderI[i]] = true;
	}

	while (!Q.empty())
	{
		int currentSurface = Q.front(); Q.pop();

		int neighborSize = dualC[currentSurface].size();
		for (int i = 0; i < neighborSize; i++)
		{
			if (dualC[currentSurface][i].priority > 0)
			{	//只要有通路
				int newSurface = dualC[currentSurface][i].order;
				if (ifIn[newSurface] == false)
				{	//就把能走到的模糊区归类为I区
					Q.push(newSurface);
					ifIn[newSurface] = true;
					block[newSurface] = BLOCK(REPI, 0, true);
				}
			}
		}
	}

	for (int i = 0; i < fnum; i++)
	{
		if (dualC[i].size() > 0)
		{	//剩下的模糊区自然是J区
			if (block[i].ifConfirm == false)
			{
				block[i].ifConfirm = true;
				block[i].area1 = REPJ;
			}
		}
	}
}






void createFlowNetwork(vector<vector<surfacePriority>>& angdDual, vector<BLOCK>& block, vector<int>& borderI, double avgAngd,
	int REPI, int REPJ)
{
	int fnum = block.size();
	vector<bool> ifDelete(fnum, false); //储存每个节点在这一轮判断中是否被整个删掉

	for (int i = 0; i < fnum; i++)
	{	//对于每个节点：
		if (block[i].ifConfirm == false)
		{	//如果是模糊区
			if ((block[i].area1 == REPI && block[i].area2 == REPJ) || (block[i].area1 == REPJ && block[i].area2 == REPI))
			{	//如果是ij模糊区,不用做任何操作
				continue;
			}
			else
			{	//如果不是ij模糊区，直接删掉
				angdDual[i].clear();
				ifDelete[i] = true;
			}
		}
		else
		{
			if (block[i].area1 == REPI || block[i].area1 == REPJ)
			{	//如果是ij区
				int neighborSize = angdDual[i].size();
				bool ifBorder = false;
				for (int j = 0; j < neighborSize; j++)
				{
					int thisOrder = angdDual[i][j].order;
					if (block[thisOrder].ifConfirm == false &&
						((block[thisOrder].area1 == REPI && block[thisOrder].area2 == REPJ) ||
							(block[thisOrder].area1 == REPJ && block[thisOrder].area2 == REPI)))
					{	//如果这个ij区面片的某个邻域是ij模糊区，那么把它放进border中
						ifBorder = true;
						break;
					}
				}

				if (!ifBorder)
				{	//如果它不是ij区与ij模糊区边界，直接删掉
					angdDual[i].clear();
					ifDelete[i] = true;
				}
				else
				{	//如果它是边界且是i的边界，存进border，如果是j的边界，什么也不用作
					if (block[i].area1 == REPI)
						borderI.push_back(i);
				}

			}
			else
			{	//如果是ij以外的其他确定区,直接删掉
				angdDual[i].clear();
				ifDelete[i] = true;
			}
		}
	}

	//将上述被删掉的边的入边也删掉
	for (int i = 0; i < fnum; i++)
	{
		if (angdDual[i].size() != 0)
		{//对于没有被删掉的边
			vector<surfacePriority>::iterator it = angdDual[i].begin();
			while (it != angdDual[i].end())
			{
				if (ifDelete[it->order])	//如果某条边指向了已经被删掉的面片
					angdDual[i].erase(it);
				else
					it++;
			}
		}
	}

	//然后将边权值由角距离改成最大水流量
	for (int i = 0; i < fnum; i++)
	{
		int size = angdDual[i].size();
		for (int j = 0; j < size; j++)
		{
			angdDual[i][j].priority = 1 / ((angdDual[i][j].priority / avgAngd) + 1);
		}
	}
}



