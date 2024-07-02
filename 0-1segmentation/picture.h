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
#define DEBUG_PRECISEDIVIDE
using namespace std;



void createDual(vector<surface>& surfaceData, vector<vector<surfacePriority> >& dual, double& avgAngd, vector<vector<surfacePriority> >& angdDual);
//dual应当是一个邻接表，因为每个surface只有三个邻边，非常稀疏
int dijkstra(vector<vector<double> >& dual, int s, int fnum, vector<double>& priority);
int dijkstra_prior(vector<vector<surfacePriority> >& dual, int s, int fnum, vector<double>& minDist);	//用优先级队列
void posibility(vector<double>& PA, vector<double>& PB, vector<vector<double> >& dist, int& REPA, int& REPB, int fnum, int times);
void modify(vector<double>& PA, vector<double>& PB, vector<vector<double> >& dist, int& REPA, int& REPB, int fnum, int times);
void findMinCut(vector<vector<surfacePriority> >& dualC, vector<int>& borderA, vector<area>& Area);
void preciseDivide(vector<vector<surfacePriority> >& dualC, vector<int>& borderA, vector<area>& Area);
void createFlowNetwork(vector<vector<surfacePriority> >& angdDual, vector<area>& Area, vector<int>& borderA, double avgAngd);

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
			if (dual[s][j]&&(!ifVisited[j]))
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

#ifndef DEBUG_PICTURE
		cout << "thisSurface: " << thisOrder << ", " << thisPriority << endl;
		cout << "now the minDist is: " << endl;
		for (int i = 0; i < fnum; i++)
			cout << minDist[i] << "  ";
		cout << endl;
		cout << "neibours: " << endl;
#endif // !DEBUG_PICTURE



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

#ifndef DEBUG_PICTURE
				cout << i << ": " << newOrder << ", " << newPriority << endl;
#endif // !DEBUG_PICTURE
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








void posibility(vector<double>& PA, vector<double>& PB, vector<vector<double> >& dist, int& REPA, int& REPB, int fnum, int times)
//其中priority是上一步中用dijkstra算出的所有i到j的最短路
{
	if (times == ITERATION_TIMES)
		return;
#ifndef DEBUG1
	cout << "posibility-turn" << times << ": REPA=" << REPA << ", REPB=" << REPB << endl;
#endif // !DEBUG1

	for (int i = 0; i < fnum; i++)
	{
	    if(i == REPA)
            PA[i] = 1;
        else
            PA[i] = (1/dist[i][REPA]) / (1/dist[i][REPA] + 1/dist[i][REPB]);

        if(i==REPB)
            PB[i] = 1;
        else
            PB[i] = (1/dist[i][REPB]) / (1/dist[i][REPA] + 1/dist[i][REPB]);
	}

    #ifndef DEBUG1
	cout<<"PA: "<<endl;
	for(int i=0;i<fnum;i++)
        cout<<PA[i]<<" ";
    cout<<endl;
    cout<<"PB: "<<endl;
	for(int i=0;i<fnum;i++)
        cout<<PB[i]<<" ";
    cout<<endl;
	#endif // DEBUG1

	modify(PA, PB, dist, REPA, REPB, fnum, times);

}

void modify(vector<double>& PA, vector<double>& PB, vector<vector<double> >& dist, int& REPA, int& REPB, int fnum, int times)
{
#ifndef DEBUG1
	cout << "modify-turn" << times << endl;
#endif // !DEBUG1

	int newREPA = REPA;
	double minSum = INT_MAX;
	for (int i = 0; i < fnum; i++)
	{
		double sumPD = 0;
		for (int j = 0; j < fnum; j++)
			sumPD += PA[j] * dist[j][i];

		if (sumPD < minSum)
		{
			minSum = sumPD;
			newREPA = i;
		}
	}

	int newREPB = REPB;
	minSum = INT_MAX;
	for (int i = 0; i < fnum; i++)
	{
		double sumPD = 0;
		for (int j = 0; j < fnum; j++)
			sumPD += PB[j] * dist[j][i];

		if (sumPD < minSum)
		{
			minSum = sumPD;
			newREPB = i;
		}
	}
    cout<<endl; cout<<endl;


	if ((newREPA == REPA && newREPB == REPB)||(newREPA == REPB && newREPB == REPA))
		return;
	else
	{
		REPA = newREPA;
		REPB = newREPB;
		posibility(PA, PB, dist, REPA, REPB, fnum, times + 1);
	}
}


void findMinCut(vector<vector<surfacePriority> >& dualC, vector<int>& borderA, vector<area>& Area)
{
    #ifndef DEBUG_FINDMINCUT
    cout<<"a new turn"<<endl;
    #endif // DEBUG_FINDMINCUT

	queue<int> Q;

	int sizeA = borderA.size();	//起始点的数量
	int fnum = Area.size();
	vector<int> pre(fnum, -1);	//用于储存每个面片在路径上的前一个面片,不在路径中的点置为-1
	vector<bool> ifIn(fnum, false);

	for (int i = 0; i < sizeA; i++)
	{
		Q.push(borderA[i]);	//先将起始点全部入队
		ifIn[borderA[i]] = true;
	}

	#ifndef DEBUG_FINDMINCUT
	cout << "border has been already in" << endl;
	#endif // DEBUG_FINDMINCUT

	bool ifFound = false;
	int end;
	while (!Q.empty() && !ifFound) //只要Q空了或者找到了B区点就结束循环
	{
		int thisSurface = Q.front(); Q.pop();
		int neighborSize = dualC[thisSurface].size();
		for (int i = 0; i < neighborSize; i++)
		{
			if (dualC[thisSurface][i].priority > 0)
			{//如果当前节点到它的第i个邻域的边权值是正的
				int newSurface = dualC[thisSurface][i].order;
				if(ifIn[newSurface] == false)
                {//并且这个邻域的节点还没有入过队
                    pre[newSurface] = thisSurface;
                    ifIn[newSurface] = true;
                    if (Area[newSurface] == AREA_B)	//如果有对面边界的面片入队
                    {
                        ifFound = true;
                        end = newSurface;
                        break;
                    }
                    Q.push(newSurface);
                }
			}
		}
	}

	if (ifFound)	//如果找到了，沿着这条路更新边权值
	{
		double minCap = INT_MAX;
		int currentSurface = end;
		while (1)
		{
			int preSurface = pre[currentSurface];

			#ifndef DEBUG_FINDMINCUT
			cout<<"->"<<preSurface;
			#endif // DEBUG_FINDMINCUT

			if (preSurface < 0)
				break;	//如果前一个面片编号为-1，证明现在这个就是起始点
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
		cout<<endl;
		cout<<"the min cap is: "<<minCap<<endl;
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

		//更新完邻接表dualC之后，继续寻找路径，直到找不到通路
		findMinCut(dualC, borderA, Area);
	}
	else
    {
        #ifndef DEBUG_FINDMINCUT
        cout<<endl;
        cout<<"terminated"<<endl;
        #endif // DEBUG_FINDMINCUT
        return;
    }

	//如果没找到通路就直接退出
}

void preciseDivide(vector<vector<surfacePriority> >& dualC, vector<int>& borderA, vector<area>& Area)
{
	queue<int> Q;

	int fnum = Area.size();
	vector<bool> ifIn(fnum, false);

	int sizeA = borderA.size();
	for (int i = 0; i < sizeA; i++)
	{	//所有初始节点入队
		Q.push(borderA[i]);
		ifIn[borderA[i]] = true;
	}

	#ifndef DEBUG_PRECISEDIVIDE
	cout<<"border has been in: "<<endl;
	int borderSize = borderA.size();
	for(int i=0;i<borderSize;i++)
    {
        cout<<borderA[i]<<"  ";
    }
    cout<<endl;
	#endif // DEBUG_PRECISEDIVIDE

	while (!Q.empty())	//这次不是找通路了，而是广度优先遍历
	{
		int currentSurface = Q.front();
		Q.pop();

		#ifndef DEBUG_PRECISEDIVIDE
		cout<<"currentDurface: "<<currentSurface<<endl;
		#endif // DEBUG_PRECIDEDIVIDE

		int neighborSize = dualC[currentSurface].size();
		for (int i = 0; i < neighborSize; i++)
		{
			if (dualC[currentSurface][i].priority > 0)
			{	//把能够走通到达的所有节点都归为A区
				int newSurface = dualC[currentSurface][i].order;
				if(ifIn[newSurface]==false)
                {
                    #ifndef DEBUG_PRECISEDIVIDE
                    cout<<"surface"<<newSurface<<" has been in"<<endl;
                    #endif // DEBUG_PRECISEDIVIDE
                    Q.push(newSurface);
                    ifIn[newSurface] = true;
                    Area[newSurface] = AREA_A;
                }
			}
		}
	}//首先将所有从A侧能遍历到的所有面片都标记为A

	for (int i = 0; i < fnum; i++)
	{
		if (Area[i] == AREA_C)
			Area[i] = AREA_B;
	}//剩下的C区自然都是B的
}

void createFlowNetwork(vector<vector<surfacePriority> >& angdDual, vector<area>& Area, vector<int>& borderA, double avgAngd)
{//根据权值为角距离的对偶图以及ABC区的划分，生成一个只包括C区和与之相邻的AB区面片的邻接表
	int fnum = Area.size();
	for (int i = 0; i < fnum; i++)
	{
		if (Area[i] == AREA_A || Area[i] == AREA_B)
		{
			int size = angdDual[i].size();
			bool ifBorder = false;
			for (int j = 0; j < size; j++)
			{
				int thisOrder = angdDual[i][j].order;
				if (Area[thisOrder] == AREA_C)
				{	//如果一个A/B区的面片的邻域中有C区面片
					ifBorder = true;
					break;
				}
			}

			if (!ifBorder)		//如果这个A,B区的面片并不与C区相邻,则断开它的所有联系
				angdDual[i].clear();
			else
			{//如果这个面片与C区相邻,
				//存入对应的边界(只用存A的)
				if (Area[i] == AREA_A)
					borderA.push_back(i);

				//然后看它的邻域，删掉所有不是C区的相邻面片
				vector<surfacePriority>::iterator it = angdDual[i].begin();
				while (it != angdDual[i].end())
				{
					if (Area[it->order] != AREA_C)
						angdDual[i].erase(it);
					else
						it++;
				}
			}
		}
	}

	//此时angdDual中该删掉的都已经删掉，剩下的边权值为两面片的角距离，需要更换成定义的最大水流量、
	for (int i = 0; i < fnum; i++)
	{
		int size = angdDual[i].size();
		for (int j = 0; j < size; j++)
		{
			angdDual[i][j].priority = 1 / ((angdDual[i][j].priority / avgAngd) + 1);
		}
	}
}

