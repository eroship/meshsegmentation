#include <iostream>
#include <set>
#include "geometry.h"
#include "picture.h"
#define EPSILON 0.1
#define MAX_CATICORY 12
#define DEBUG
using namespace std;
int main()
{

	string inFile;
	cin >> inFile;
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

	//构建邻接表
	vector<vector<surfacePriority> > dual(fnum);
	vector<vector<surfacePriority> > angdDual(fnum);
	double avgAngd = 0;
	createDual(surfaceData, dual, avgAngd, angdDual);	
	cout << "dual created" << endl;
	//除了标准的邻接表之外，还保存一下角距离为权值的邻接表和平均角距离，后面有用

	//第二步：生成聚类种子
	//思路：还是先算出两两之间的最短距离
	//先找距离最远的两个点，再找到：到它们两个点的距离之和最远的点
	//并记录到它们两个的距离的最小值（G）
	//以此类推再往上加点，找到G下降最快的时候
	vector<vector<double> > priority(fnum, vector<double>(fnum, 0));	//储存两两之间的最短距离
	for (int i = 0; i < fnum; i++)
	{
		dijkstra_prior(dual, i, fnum, priority[i]);
	}

	
	vector<int> REP;	//储存聚类中心，即所有是聚类中心的面片的下标
	vector<double> REPMinDist;	//新的聚类中心到之前的所有聚类中心的距离的最小值
	//此向量的第i个表示有i个聚类中心的距离最小值
	//先做两个的，即找两个距离最远的
	double maxDist = 0;
	int REPA = 0, REPB = 0;
	for (int i = 0; i < fnum; i++)
	{
		for (int j = 0; j < fnum; j++)
		{
			if (priority[i][j] > maxDist)
			{
				maxDist = priority[i][j];
				REPA = i;
				REPB = j;
			}
		}
	}
	REP.push_back(REPA); REP.push_back(REPB); 
	REPMinDist.push_back(0); REPMinDist.push_back(0);	//当有0个或者一个聚类中心时，距离都是0
	REPMinDist.push_back(maxDist);	//第二个聚类中心到第一个聚类中心的距离
	cout << "the first two REP: " << REPA << " " << REPB << endl;
	
	int knum = 3;	//从三个聚类中心开始算起
	for (int i = 3; i <= MAX_CATICORY; i++)
	{
		int newREP = 0;
		double maxDist = 0;
		for (int k = 0; k < fnum; k++)
		{//对于第k个点
			double sumDist = 0;
			for (int j = 0; j < knum - 1; j++)
				sumDist += priority[k][REP[j]]; //第k个点到第j个REP的距离

			if (sumDist > maxDist)
			{	//找到离前面的所有中心距离之和最大的点
				maxDist = sumDist;
				newREP = k;
			}
		}
		cout << "REP" << i << ": " << newREP << endl;

		//得到newREP之后，计算它与前面的REP的最小距离
		double minDist = INT_MAX;
		for (int j = 0; j < knum - 1; j++)
		{
			cout << priority[newREP][REP[j]] << "  ";
			if (priority[newREP][REP[j]] < minDist)
				minDist = priority[newREP][REP[j]];
		}
		cout << "minDist: " << minDist << endl;

		//新的聚类中心和最小距离都得到了，将它们放进对应的向量中
		REP.push_back(newREP);
		REPMinDist.push_back(minDist);
		knum++;
	}


	cout << "the G function of each number of REP: " << endl;
	for (int i = 0; i < knum; i++)
	{
		cout << REPMinDist[i] << "  " << endl;
	}
	cout << endl;

	//看REPMinDist中哪两个元素下降最快
	double maxDMinDist = 0;
	for (int i = 0; i < REPMinDist.size()-1; i++)
	{
		double dMinDist = REPMinDist[i] - REPMinDist[i + 1];
		if (dMinDist > maxDMinDist)
		{
			maxDMinDist = dMinDist;
			knum = i;
		}
	}
	REP.resize(knum);
	
#ifndef DEBUG
	cout << "knum: " << knum << endl;
	cout << "REP: ";
	for (int i = 0; i < knum; i++)
		cout << REP[i] << " " << endl;
	cout << endl;
#endif // !DEBUG


	
	
	//第三步：聚类
	//与2路分解一样，迭代进行概率计算和聚类中心的更新
	vector<vector<double>> psb(knum, vector<double>(fnum, 0));
	//概率矩阵：knum个向量，每个向量表示对应的点属于对应的类的概率
	posibility(psb, priority, REP, knum, fnum, 0);



	//第四步：模糊分割：
	//将某个聚类中心概率大于0.5+epsilon的面片归类到此中心
	//将不确定的面片按照概率最大的两个概率中心划分

	vector<BLOCK> block(fnum); //储存每个面片的区域信息（BLOCK型）
	for (int i = 0; i < fnum; i++)
	{
		int maxREP = 0; double maxPSB = 0;
		//先看第i个面片概率最大的聚类中心的概率是否超过了0.5+epsilon
		for (int j = 0; j < knum; j++)
		{
			if (psb[j][i] > maxPSB)
			{
				maxPSB = psb[j][i];
				maxREP = j;
			}
		} //注意概率最大的聚类中心在REP中的下标是maxREP，它在surface中的下标不是这个

		if (maxPSB > 0.5 + EPSILON)
			block[i] = BLOCK(REP[maxREP], 0, true); //确定的话，分区的第二中心默认为0
		else
		{	//否则找概率第二大的聚类中心
			int subMaxREP = 0; double subMaxPSB = 0;
			for (int j = 0; j < knum; j++)
			{
				if (j == maxREP)
					continue;	//遇到概率最大的那个聚类中心就直接跳过
				else
				{
					if (psb[j][i] > subMaxPSB)
					{
						subMaxPSB = psb[j][i];
						subMaxREP = j;
					}
				}
			}

			block[i] = BLOCK(REP[maxREP], REP[subMaxREP], false);
		}
	}
	vector<int> colorHash(fnum, -1);
	for (int i = 0; i < knum; i++)
	{
		colorHash[REP[i]] = i;	//直接访问聚类中心得到下标
	}
	outputObj(outFile1, block, pointData, surfaceData, knum, colorHash);
	
	//第五步：精细分割
	for (int i = 0; i < knum; i++)
	{
		for (int j = 0; j < i; j++)
		{//每一步，只留下ij的模糊区以及模糊区与i，j交界的区域，在此循环中j始终小于i  
			vector<int> borderI;
			vector<vector<surfacePriority>> angdDualIJ = angdDual;
			//生成一个只包含ij边缘以及模糊区的流网络angdDualIJ
			createFlowNetwork(angdDualIJ, block, borderI, avgAngd, REP[i], REP[j]);


			//找到angdDualIJ的最小割
			findMinCut(angdDualIJ, borderI, REP[j], block);

			//精细划分
			preciseDivide(angdDualIJ, borderI, REP[i], REP[j], block);
		}
	}

	outputObj(outFile2, block, pointData, surfaceData, knum, colorHash);

}


