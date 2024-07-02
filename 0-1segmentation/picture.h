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
//dualӦ����һ���ڽӱ���Ϊÿ��surfaceֻ�������ڱߣ��ǳ�ϡ��
int dijkstra(vector<vector<double> >& dual, int s, int fnum, vector<double>& priority);
int dijkstra_prior(vector<vector<surfacePriority> >& dual, int s, int fnum, vector<double>& minDist);	//�����ȼ�����
void posibility(vector<double>& PA, vector<double>& PB, vector<vector<double> >& dist, int& REPA, int& REPB, int fnum, int times);
void modify(vector<double>& PA, vector<double>& PB, vector<vector<double> >& dist, int& REPA, int& REPB, int fnum, int times);
void findMinCut(vector<vector<surfacePriority> >& dualC, vector<int>& borderA, vector<area>& Area);
void preciseDivide(vector<vector<surfacePriority> >& dualC, vector<int>& borderA, vector<area>& Area);
void createFlowNetwork(vector<vector<surfacePriority> >& angdDual, vector<area>& Area, vector<int>& borderA, double avgAngd);

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

#ifndef DEBUG_PICTURE
		cout << "thisSurface: " << thisOrder << ", " << thisPriority << endl;
		cout << "now the minDist is: " << endl;
		for (int i = 0; i < fnum; i++)
			cout << minDist[i] << "  ";
		cout << endl;
		cout << "neibours: " << endl;
#endif // !DEBUG_PICTURE



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
//����priority����һ������dijkstra���������i��j�����·
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

	int sizeA = borderA.size();	//��ʼ�������
	int fnum = Area.size();
	vector<int> pre(fnum, -1);	//���ڴ���ÿ����Ƭ��·���ϵ�ǰһ����Ƭ,����·���еĵ���Ϊ-1
	vector<bool> ifIn(fnum, false);

	for (int i = 0; i < sizeA; i++)
	{
		Q.push(borderA[i]);	//�Ƚ���ʼ��ȫ�����
		ifIn[borderA[i]] = true;
	}

	#ifndef DEBUG_FINDMINCUT
	cout << "border has been already in" << endl;
	#endif // DEBUG_FINDMINCUT

	bool ifFound = false;
	int end;
	while (!Q.empty() && !ifFound) //ֻҪQ���˻����ҵ���B����ͽ���ѭ��
	{
		int thisSurface = Q.front(); Q.pop();
		int neighborSize = dualC[thisSurface].size();
		for (int i = 0; i < neighborSize; i++)
		{
			if (dualC[thisSurface][i].priority > 0)
			{//�����ǰ�ڵ㵽���ĵ�i������ı�Ȩֵ������
				int newSurface = dualC[thisSurface][i].order;
				if(ifIn[newSurface] == false)
                {//�����������Ľڵ㻹û�������
                    pre[newSurface] = thisSurface;
                    ifIn[newSurface] = true;
                    if (Area[newSurface] == AREA_B)	//����ж���߽����Ƭ���
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

	if (ifFound)	//����ҵ��ˣ���������·���±�Ȩֵ
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
				break;	//���ǰһ����Ƭ���Ϊ-1��֤���������������ʼ��
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
		cout<<endl;
		cout<<"the min cap is: "<<minCap<<endl;
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

		//�������ڽӱ�dualC֮�󣬼���Ѱ��·����ֱ���Ҳ���ͨ·
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

	//���û�ҵ�ͨ·��ֱ���˳�
}

void preciseDivide(vector<vector<surfacePriority> >& dualC, vector<int>& borderA, vector<area>& Area)
{
	queue<int> Q;

	int fnum = Area.size();
	vector<bool> ifIn(fnum, false);

	int sizeA = borderA.size();
	for (int i = 0; i < sizeA; i++)
	{	//���г�ʼ�ڵ����
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

	while (!Q.empty())	//��β�����ͨ·�ˣ����ǹ�����ȱ���
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
			{	//���ܹ���ͨ��������нڵ㶼��ΪA��
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
	}//���Ƚ����д�A���ܱ�������������Ƭ�����ΪA

	for (int i = 0; i < fnum; i++)
	{
		if (Area[i] == AREA_C)
			Area[i] = AREA_B;
	}//ʣ�µ�C����Ȼ����B��
}

void createFlowNetwork(vector<vector<surfacePriority> >& angdDual, vector<area>& Area, vector<int>& borderA, double avgAngd)
{//����ȨֵΪ�Ǿ���Ķ�żͼ�Լ�ABC���Ļ��֣�����һ��ֻ����C������֮���ڵ�AB����Ƭ���ڽӱ�
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
				{	//���һ��A/B������Ƭ����������C����Ƭ
					ifBorder = true;
					break;
				}
			}

			if (!ifBorder)		//������A,B������Ƭ������C������,��Ͽ�����������ϵ
				angdDual[i].clear();
			else
			{//��������Ƭ��C������,
				//�����Ӧ�ı߽�(ֻ�ô�A��)
				if (Area[i] == AREA_A)
					borderA.push_back(i);

				//Ȼ����������ɾ�����в���C����������Ƭ
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

	//��ʱangdDual�и�ɾ���Ķ��Ѿ�ɾ����ʣ�µı�ȨֵΪ����Ƭ�ĽǾ��룬��Ҫ�����ɶ�������ˮ������
	for (int i = 0; i < fnum; i++)
	{
		int size = angdDual[i].size();
		for (int j = 0; j < size; j++)
		{
			angdDual[i][j].priority = 1 / ((angdDual[i][j].priority / avgAngd) + 1);
		}
	}
}

