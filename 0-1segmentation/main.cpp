#include <iostream>
#include "geometry.h"
#include "picture.h"
#define EPSILON 0.1
#define DEBUG
using namespace std;
int main()
{

	string inFile;
	cin>>inFile;
	string outFile1 = "output1.obj";
	string outFile2 = "output2.obj";

	//第一步：构建对偶图（dual）
	//对偶图用邻接矩阵表示,每条边的权重是两个节点（面片）的距离

	vector<point> pointData;
	vector<surface> surfaceData;

	getData(inFile, pointData, surfaceData);	//先获取数据
	cout << "data confirm" << endl;
	int vnum = pointData.size();
	int fnum = surfaceData.size();

	vector<vector<surfacePriority> > dual(fnum);	//这是一个邻接表
	vector<vector<surfacePriority> > angdDual(fnum);	//邻接表，储存角距离参数
	double avgAngd = 0;

	createDual(surfaceData, dual, avgAngd, angdDual);
	cout << "dual created" << endl;

	//第二步：生成聚类种子
	int start = 0, end = 0;
	double maxDist = 0;
	vector<vector<double> > priority(fnum, vector<double>(fnum, 0));
	for (int i = 0; i < fnum; i++)
	{//用priority储存两两之间的最短路，每次更新一行
		int thisEnd = dijkstra_prior(dual, i, fnum, priority[i]);
		if (priority[i][thisEnd] > maxDist)
		{
			maxDist = priority[i][thisEnd];
			start = i;
			end = thisEnd;
		}
	}

	//测试：
	cout << endl;
	//cout << "priority matrix: " << endl;
	//for (int i = 0; i < fnum; i++)
	//{
	//	for (int j = 0; j < fnum; j++)
	//		cout << setprecision(3) << setw(6) << priority[i][j];
	//	cout << endl;
	//}
	cout << "the farest pair: " << start << ", " << end << endl;

	//聚类：
	int REPA = start;
	int REPB = end;

	vector<double> PA(fnum, 0);
	vector<double> PB(fnum, 0);
	posibility(PA, PB, priority, REPA, REPB, fnum, 0);

	//模糊分割
	vector<int> CA;
	vector<int> CB;
	vector<int> CC;	//这三个向量储存三个区的面片的索引下标
	vector<area> Area(fnum);//这个向量储存每个面片位于哪个区
	for (int i = 0; i < fnum; i++)
	{
		if (PA[i] > EPSILON + 0.5)
		{
			CA.push_back(i);
			Area[i] = AREA_A;
		}

		else
		{
			if (PB[i] > EPSILON + 0.5)
			{
				CB.push_back(i);
				Area[i] = AREA_B;
			}
			else
			{
				CC.push_back(i);
				Area[i] = AREA_C;
			}
		}
	}


    outputObj(outFile1, Area, pointData, surfaceData);
#ifndef DEBUG
	outputObj(outFile, Area, pointData, surfaceData);
#endif // !DEBUG


	//精细分割：最大流
	vector<int> borderA(0);	//储存那些面片是AC相邻的
	createFlowNetwork(angdDual, Area, borderA, avgAngd);
	cout << "flow network created" << endl;

	#ifndef DEBUG
	cout<<"matrix angdDual: "<<endl;
	for(int i=0;i<fnum;i++)
    {
        int size = angdDual[i].size();
        cout<<"surface"<<i<<": ";
        for(int j=0;j<size;j++)
        {
            cout<<angdDual[i][j].order<<"("<<angdDual[i][j].priority<<")  ";
        }
        cout<<endl;
    }
	#endif // DEBUG

	//根据新获得的邻接表找最大流/最小割,得到的是一个已经出现最小割集的邻接表
	findMinCut(angdDual, borderA, Area);
	cout << "min cut created" << endl;

	#ifndef DEBUG
	cout<<"matrix angdDual: "<<endl;
	for(int i=0;i<fnum;i++)
    {
        int size = angdDual[i].size();
        cout<<"surface"<<i<<": ";
        for(int j=0;j<size;j++)
        {
            cout<<angdDual[i][j].order<<"("<<angdDual[i][j].priority<<")  ";
        }
        cout<<endl;
    }
	#endif // DEBUG

	//利用这个邻接表将C区精确划分
	preciseDivide(angdDual, borderA, Area);
	cout << "division created" << endl;

    outputObj(outFile2, Area, pointData, surfaceData);
}

