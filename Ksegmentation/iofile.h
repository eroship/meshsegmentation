#pragma once
//头文件ioFile，用于文件读取和输出
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "picture.h"
#include "geometry.h"
#define DEBUG_IOFILE
struct BLOCK
{
	int area1;
	int area2; //概率最大的两个聚类中心作为面片的下标
	bool ifConfirm;
	BLOCK(int a1 = 0, int a2 = 0, bool ifCf = false)
	{	//始终保证area1的值大于等于area2的值
		if (a1 > a2)
		{
			area1 = a1;
			area2 = a2;
		}
		else
		{
			area1 = a2;
			area2 = a1;
		}
		ifConfirm = ifCf;
	}
};
//这个结构体用于保存某个面片属于的区域信息，area1和area2表示两个概率最大的聚类中心

using namespace std;

const string COLOR[12] = { 
"255 0 0", "0 255 0", "0 0 255", "124 124 0", "124 0 124", "0 124 124",
"0 124 255", "0 255 124", "124 0 255", "124 255 0", "255 124 0", "255 0 124"};

void getData(string fileName, vector<point>& pointData, vector<surface>& surfaceData);
double stringToDouble(string s);
//void outputObj(string filename, vector<area>& Area, vector<point>& pointData, vector<surface>& surfaceData);
void outputObj(string filename, vector<BLOCK>& block, vector<point>& pointDta, vector<surface>& surfaceData, int knum, vector<int>& colorHash);


void getData(string fileName, vector<point>& pointData, vector<surface>& surfaceData)
{
	//文件输入流，打开文件
	ifstream streamIn(fileName.data());
	if (!streamIn.is_open())
	{
		cout << "open failed" << endl;
		exit(0);
	}

	string s;
	string s0, s1, s2, s3, s4;

	while (getline(streamIn, s))
	{
		if (s[0] == '#')
			continue;

		//字符串输入流，读取到的数据存到s1~中
		istringstream stringData(s);
		stringData >> s0 >> s1 >> s2 >> s3;
		if (s[0] == 'v')
		{
			//pointData记录所有的点
			point thisPoint(stringToDouble(s1), stringToDouble(s2), stringToDouble(s3));
			pointData.push_back(thisPoint);

#ifndef DEBUG_IOFILE
			cout << "v:" << thisPoint.x << ", " << thisPoint.y << ", " << thisPoint.z << endl;
			cout << endl;
#endif // !DEBUG

		}
		if (s[0] == 'f')
		{
			int pt1 = (int)stringToDouble(s1);
			int pt2 = (int)stringToDouble(s2);
			int pt3 = (int)stringToDouble(s3);
			surface thisSurface(pointData[pt1 - 1], pointData[pt2 - 1], pointData[pt3 - 1], pt1 - 1, pt2 - 1, pt3 - 1);
			surfaceData.push_back(thisSurface);

#ifndef DEBUG_IOFILE
			cout << "f: " << thisSurface.vtx1 << ", " << thisSurface.vtx2 << ", " << thisSurface.vtx3 << endl;
			cout << "center: (" << thisSurface.centerx << ", " << thisSurface.centery << ", " << thisSurface.centerz << ")" << endl;
			cout << "normal vector: (" << thisSurface.nmvx << ", " << thisSurface.nmvy << ", " << thisSurface.nmvz << ")" << endl;
			cout << endl;
#endif // !DEBUG

		}
	}
}

double stringToDouble(string s)
{
	int mul = 1;
	double sum = 0;
	int rp = 1;
	for (int t = s.size() - 1; t >= 0; --t)
	{
		if (s[t] == '.')
			rp = mul;
		else
		{
			if (s[t] == '-') rp = -rp;
			else
			{
				sum += (double)(mul * ((int)s[t] - 48));
				mul *= 10;
			}
		}
	}
	return sum / rp;
}



void outputObj(string filename, vector<BLOCK>& block, vector<point>& pointData, vector<surface>& surfaceData, int knum, vector<int>& colorHash)
{
	int vnum = pointData.size();
	int fnum = surfaceData.size();
	vector<int> pointColor(vnum);

	for (int i = 0; i < fnum; i++)
	{
		//如果是未确定区，不涂颜色
		if (block[i].ifConfirm == false)
		{
			pointColor[surfaceData[i].vtx1] = -1;
			pointColor[surfaceData[i].vtx2] = -1;
			pointColor[surfaceData[i].vtx3] = -1;
		}
		else
		{
			int color = colorHash[block[i].area1];
			pointColor[surfaceData[i].vtx1] = color;
			pointColor[surfaceData[i].vtx2] = color;
			pointColor[surfaceData[i].vtx3] = color;
		}

	}

	ofstream output(filename.data());
	for (int i = 0; i < vnum; i++)
	{
		output << "v " << pointData[i].x << " " << pointData[i].y << " " << pointData[i].z << " ";
		if (pointColor[i] == -1)
			output << "255 255 255" << endl;
		else
			output << COLOR[pointColor[i]] << endl;
	}

	for (int i = 0; i < fnum; i++)
	{
		output << "f " << surfaceData[i].vtx1 + 1 << " " << surfaceData[i].vtx2 + 1 << " " << surfaceData[i].vtx3 + 1 << endl;
	}
}

