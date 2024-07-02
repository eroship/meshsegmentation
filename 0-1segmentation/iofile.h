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
typedef enum { AREA_A, AREA_B, AREA_C } area;

using namespace std;

#define COLORA "255 0 0"
#define COLORB "0 255 0"
#define COLORC "0 0 255"

void getData(string fileName, vector<point>& pointData, vector<surface>& surfaceData);
double stringToDouble(string s);
void outputObj(string filename, vector<area>& Area, vector<point>& pointData, vector<surface>& surfaceData);

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

void outputObj(string filename, vector<area>& Area, vector<point>& pointData, vector<surface>& surfaceData)
{
	int vnum = pointData.size();
	int fnum = surfaceData.size();
	vector<area> color(vnum);	//储存每个点位于什么区，交界处的点随意

	for (int i = 0; i < fnum; i++)
	{
		switch (Area[i])
		{
		case AREA_A:
			color[surfaceData[i].vtx1] = AREA_A;
			color[surfaceData[i].vtx2] = AREA_A;
			color[surfaceData[i].vtx3] = AREA_A;
			break;
		case AREA_B:
			color[surfaceData[i].vtx1] = AREA_B;
			color[surfaceData[i].vtx2] = AREA_B;
			color[surfaceData[i].vtx3] = AREA_B;
			break;
		case AREA_C:
			color[surfaceData[i].vtx1] = AREA_C;
			color[surfaceData[i].vtx2] = AREA_C;
			color[surfaceData[i].vtx3] = AREA_C;
			break;
		default:
			break;
		}
	}

	ofstream output(filename.data());
	for (int i = 0; i < vnum; i++)
	{
		output << "v " << pointData[i].x << " " << pointData[i].y << " " << pointData[i].z << " ";
		switch (color[i])
		{
		case AREA_A: output << COLORA << endl; break;
		case AREA_B: output << COLORB << endl; break;
		case AREA_C: output << COLORC << endl; break;
		default:
			break;
		}
	}

	for (int i = 0; i < fnum; i++)
	{
		output << "f " << surfaceData[i].vtx1 + 1 << " " << surfaceData[i].vtx2 + 1 << " " << surfaceData[i].vtx3 + 1 << endl;
	}

}
