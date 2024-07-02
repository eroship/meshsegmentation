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
int dijkstra_prior(vector<vector<surfacePriority> >& dual, int s, int fnum, vector<double>& minDist);	//�����ȼ�����

void posibility(vector<vector<double>>& psb, vector<vector<double>>& dist, vector<int>& REP, int knum, int fnum, int times);
void modify(vector<vector<double>>& psb, vector<vector<double>>& dist, vector<int>& REP, int knum, int fnum, int times);

void findMinCut(vector<vector<surfacePriority>>& dualC, vector<int>& borderI, int REPJ, vector<BLOCK>& block);
void preciseDivide(vector<vector<surfacePriority>>& dualC, vector<int>& borderI, int REPI, int REPJ, vector<BLOCK>& block);
void createFlowNetwork(vector<vector<surfacePriority>>& angdDual, vector<BLOCK>& block, vector<int>& borderI, double avgAngd,
	int REPI, int REPJ);


void createDual(vector<surface>& surfaceData, vector<vector<surfacePriority> >& dual, double& avgAngd, vector<vector<surfacePriority> >& angdDual)
{
	map<side, int> sideHash;	//�����Ǳߣ���ֵ��ĳ����Ƭ�ı��
	int fnum = surfaceData.size();

	vector<vector<double> > dualMatrix(fnum, vector<double>(fnum, 0));	//��ʱ�����ڽӾ��󣬺�����ת���ڽӱ�

	int sideNum = 0;
	double totalGeod = 0;
	double totalAngd = 0;

	vector<vector<double> > geodMatrix(fnum, vector<double>(fnum, 0));
	vector<vector<double> > angdMatrix(fnum, vector<double>(fnum, 0));
	//����������󴢴��ؾ���ͽǾ���

	for (int i = 0; i < fnum; i++)	//��һ�ֱ������ȼ�¼�����ܺͺͱߵ�����
	{
		//��i����Ƭ��������
		side s1(surfaceData[i].vtx1, surfaceData[i].vtx2);
		side s2(surfaceData[i].vtx2, surfaceData[i].vtx3);
		side s3(surfaceData[i].vtx3, surfaceData[i].vtx1);

		pair<map<side, int>::iterator, bool> insertPair;//�������ʧ�ܣ���ô����������ļ�ֵ���Ǵ����������

#ifndef DEBUG_PICTURE
		cout << "surface" << i << ": " << endl;
		cout << "side1: " << s1.a << "-" << s1.b << ", side2: " << s2.a << "-" << s2.b << ", side3: " << s3.a << "-" << s3.b << endl;
#endif // !DEBUG_PICTURE

		insertPair = sideHash.insert(pair<side, int>(s1, i));
		if (!insertPair.second)
		{	//������ʧ�ܣ���s1����һ����ʹ��棨i�����б�s1
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
	avgAngd = totalAngd / sideNum;	//����ƽ���Ǿ���

	for (int i = 0; i < fnum; i++)
	{
		for (int j = 0; j < fnum; j++)
			dualMatrix[i][j] = theta * geodMatrix[i][j] / avgGeod + (1 - theta) * angdMatrix[i][j] / avgAngd;
	}

	//�����ڽӱ�
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
	//������ʾ�Ǿ�����ڽӱ�
	for (int i = 0; i < fnum; i++)
	{
		int size = dual[i].size();
		for (int j = 0; j < size; j++)
		{
			int thisOrder = dual[i][j].order;	//ָ�����Ƭ�ı��
			double angd = angdMatrix[i][thisOrder];//�ҵ���Ӧ�ĽǾ���
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
//�ڶ�żͼdual�У���sΪԴ�㣬Ѱ�����·��fnum�Ǿ����ά��, priority��һ�����飬����ÿ��������·�ĳ���
{
	vector<bool> ifVisited(fnum, false);	//���ڱ���Ѿ���·�����еĽڵ�

	priority[s] = 0;
	for (int i = 0; i < fnum; i++)
	{
		ifVisited[s] = true;	//ÿһ����s��Ȩֵ��С�ĵ�
		for (int j = 0; j < fnum; j++)
		{	//������s���ڵĵ㣬����û��visit�ľ͸���Ȩֵ
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
{//dualΪԭͼ���ڽӱ�priority����Դ��s�����е����̾���,���ش�s��������̾�����Զ���Ǹ���Ƭ�ı��

	priority_queue<surfacePriority> neighbor;
	minDist[s] = 0;
	vector<bool> ifVisited(fnum, false);

	neighbor.push(surfacePriority(0, s));

	while (!neighbor.empty())
	{
		surfacePriority thisSurface = neighbor.top();
		neighbor.pop();
		if (ifVisited[thisSurface.order])
			continue;	//�³��ӵĶ�ͷ�����Ѿ����ʹ��ˣ���Ϊ�ǿ����ظ���ӵ�
		ifVisited[thisSurface.order] = true;

		int thisOrder = thisSurface.order;
		double thisPriority = thisSurface.priority;
		//��ǰ���ӵĽڵ㣨����ǰȷ�Ͼ���Ľڵ㣩����Ϣ

		minDist[thisOrder] = thisPriority;	//��ÿ���ڵ����̾���浽minDist��

		//�������������򲢷������
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
{//���룺�洢���ʵľ��󣬴洢����ľ��󣬴洢�������ĵ��������������ĸ�������Ƭ�����������Ĵ���
	if (times == ITERATION_TIMES)
		return;
	for (int i = 0; i < fnum; i++)
	{//���ڵ�i���㣬���������ڵ�j�����ĵĸ���
		double sumDist = 0;
		for (int j = 0; j < knum; j++)
			sumDist += (1 / dist[i][REP[j]]); //�ȵõ������������ĵ��ܾ���

		for (int j = 0; j < knum; j++)
			psb[j][i] = (1 / dist[REP[j]][i]) / sumDist;
		//��j�����ĵ�ĸ�������
	}

	for (int i = 0; i < knum; i++)
		psb[i][REP[i]] = 1;

	modify(psb, dist, REP, knum, fnum, times);

}



void modify(vector<vector<double>>& psb, vector<vector<double>>& dist, vector<int>& REP, int knum, int fnum, int times)
{
	cout << "modify: turn" << times << endl;

	vector<int> newREP(knum);
	for (int i = 0; i < knum; i++) //һ�θ������еľ�������
	{//���ڵ�i����������
		int thisNewREP = REP[i];
		double minSum = (double)INT_MAX;
		for (int j = 0; j < fnum; j++)
		{//���ڵ�j����Ƭ��ÿ����Ƭ���м��㵽���е���Ƭ��PD������*���ʣ�֮��
			double sumPD = 0.0;
			for (int k = 0; k < fnum; k++)
				sumPD += psb[i][k] * dist[k][j];

			if (sumPD < minSum)
			{	//�ҵ���������i��Ӧ��PD����С���Ǹ��㣨�µľ������ģ�
				minSum = sumPD;
				thisNewREP = j;
			}
		}

		newREP[i] = thisNewREP;
	}

	vector<int> sREP = REP; //Ϊ���������ʱ���ƻ�ԭ�е�REP��˳���¸���һ��һ����vector
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

	vector<int> pre(fnum, -1); //����ÿ����Ƭ��·���ϵ�ǰһ���ڵ�
	vector<bool> ifIn(fnum, false);

	for (int i = 0; i < sizeI; i++)
	{	//�߽������
		Q.push(borderI[i]);
		ifIn[borderI[i]] = true;
	}

	bool ifFound = false;
	int end;
	while (!Q.empty() && !ifFound)
	{	//ֻҪ�ҵ�ͨ·���߶����ѿվͽ���ѭ��
		int thisSurface = Q.front(); Q.pop();
		int neighborSize = dualC[thisSurface].size();
		for (int i = 0; i < neighborSize; i++)
		{	//���ѳ��ӵ���Ƭ������ 
			if (dualC[thisSurface][i].priority > 0)
			{	//�����ͨ·
				int newSurface = dualC[thisSurface][i].order;
				if (ifIn[newSurface] == false)
				{	//���û�����
					pre[newSurface] = thisSurface;
					ifIn[newSurface] = true;
					if (block[newSurface].ifConfirm == true && block[newSurface].area1 == REPJ)
					{	//�������µ���Ƭ���յ�
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
	{	//����ҵ�
		double minCap = INT_MAX; //·����Сˮ����
		int currentSurface = end;
		while (1)
		{	//������·���ϵ���Сˮ����
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
				{//�ҵ�ǰ��Ƭ��ǰһ����Ƭ֮��ı�,Ȼ�����ı�Ȩֵ�Ƿ�������·����С��
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
		{	//�ҵ�����ͨ·�ϵ���С��Ȩֵ�����б�Ȩֵ��ȥ��
			int preSurface = pre[currentSurface];
			if (preSurface < 0)
				break;

			int size = dualC[preSurface].size();
			for (int i = 0; i < size; i++)
			{
				if (dualC[preSurface][i].order == currentSurface)//��ǰһ����Ƭָ��ǰ��Ƭ�ı�
				{
					dualC[preSurface][i].priority -= minCap;
					break;
				}
			}

			size = dualC[currentSurface].size();
			for (int i = 0; i < size; i++)
			{
				if (dualC[currentSurface][i].order == preSurface)//�ҵ�ǰ��Ƭָ��ǰһ����Ƭ�ı�
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
	{	//�߽������
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
			{	//ֻҪ��ͨ·
				int newSurface = dualC[currentSurface][i].order;
				if (ifIn[newSurface] == false)
				{	//�Ͱ����ߵ���ģ��������ΪI��
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
		{	//ʣ�µ�ģ������Ȼ��J��
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
	vector<bool> ifDelete(fnum, false); //����ÿ���ڵ�����һ���ж����Ƿ�����ɾ��

	for (int i = 0; i < fnum; i++)
	{	//����ÿ���ڵ㣺
		if (block[i].ifConfirm == false)
		{	//�����ģ����
			if ((block[i].area1 == REPI && block[i].area2 == REPJ) || (block[i].area1 == REPJ && block[i].area2 == REPI))
			{	//�����ijģ����,�������κβ���
				continue;
			}
			else
			{	//�������ijģ������ֱ��ɾ��
				angdDual[i].clear();
				ifDelete[i] = true;
			}
		}
		else
		{
			if (block[i].area1 == REPI || block[i].area1 == REPJ)
			{	//�����ij��
				int neighborSize = angdDual[i].size();
				bool ifBorder = false;
				for (int j = 0; j < neighborSize; j++)
				{
					int thisOrder = angdDual[i][j].order;
					if (block[thisOrder].ifConfirm == false &&
						((block[thisOrder].area1 == REPI && block[thisOrder].area2 == REPJ) ||
							(block[thisOrder].area1 == REPJ && block[thisOrder].area2 == REPI)))
					{	//������ij����Ƭ��ĳ��������ijģ��������ô�����Ž�border��
						ifBorder = true;
						break;
					}
				}

				if (!ifBorder)
				{	//���������ij����ijģ�����߽磬ֱ��ɾ��
					angdDual[i].clear();
					ifDelete[i] = true;
				}
				else
				{	//������Ǳ߽�����i�ı߽磬���border�������j�ı߽磬ʲôҲ������
					if (block[i].area1 == REPI)
						borderI.push_back(i);
				}

			}
			else
			{	//�����ij���������ȷ����,ֱ��ɾ��
				angdDual[i].clear();
				ifDelete[i] = true;
			}
		}
	}

	//��������ɾ���ıߵ����Ҳɾ��
	for (int i = 0; i < fnum; i++)
	{
		if (angdDual[i].size() != 0)
		{//����û�б�ɾ���ı�
			vector<surfacePriority>::iterator it = angdDual[i].begin();
			while (it != angdDual[i].end())
			{
				if (ifDelete[it->order])	//���ĳ����ָ�����Ѿ���ɾ������Ƭ
					angdDual[i].erase(it);
				else
					it++;
			}
		}
	}

	//Ȼ�󽫱�Ȩֵ�ɽǾ���ĳ����ˮ����
	for (int i = 0; i < fnum; i++)
	{
		int size = angdDual[i].size();
		for (int j = 0; j < size; j++)
		{
			angdDual[i][j].priority = 1 / ((angdDual[i][j].priority / avgAngd) + 1);
		}
	}
}



