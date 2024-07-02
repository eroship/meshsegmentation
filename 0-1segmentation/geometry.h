#pragma once
//头文件geometry.h
//用于一些几何和拓扑元素的定义
#include <math.h>
#define eta 0.1
#define theta 0.3
struct point
{
	double x, y, z;
	double r, g, b;

	point(double x = 0, double y = 0, double z = 0, double r = 0, double g = 0, double b = 0) :
		x(x), y(y), z(z), r(r), g(g), b(b) {}
};

struct surface
{
	point vtxp1, vtxp2, vtxp3;	//三个顶点
	int vtx1, vtx2, vtx3;//三个顶点的编号，注意是在vector中的编号而不是文件中的
	double centerx, centery, centerz;			//质心
	double nmvx, nmvy, nmvz;	//法向量
	surface(point vtxp1, point vtxp2, point vtxp3, int vtx1, int vtx2, int vtx3) :
		vtxp1(vtxp1), vtxp2(vtxp2), vtxp3(vtxp3), vtx1(vtx1), vtx2(vtx2), vtx3(vtx3)
	{
		double x12 = vtxp2.x - vtxp1.x;
		double y12 = vtxp2.y - vtxp1.y;
		double z12 = vtxp2.z - vtxp1.z;

		double x13 = vtxp3.x - vtxp1.x;
		double y13 = vtxp3.y - vtxp1.y;
		double z13 = vtxp3.z - vtxp1.z;

		nmvx = y12 * z13 - y13 * z12;
		nmvy = x12 * z13 - x13 * z12;
		nmvz = x12 * y13 - x13 * y12;

		centerx = (vtxp1.x + vtxp2.x + vtxp3.x) / 3;
		centery = (vtxp1.y + vtxp2.y + vtxp3.y) / 3;
		centerz = (vtxp1.z + vtxp2.z + vtxp3.z) / 3;
	}
};

struct surfacePriority
{
	double priority;
	int order;
	surfacePriority(double priority, int order) :priority(priority), order(order) {}
	friend bool operator<(surfacePriority sa, surfacePriority sb)
	{
		return sa.priority > sb.priority;
	}
};
//此结构用于dijkstra排序，因为使用了优先级队列，根据priority排序

struct side
{
	int a, b;

	side(int pa, int pb)
	{
		if (pa < pb)
		{
			a = pa;
			b = pb;
		}
		else
		{
			a = pb;
			b = pa;
		}
	}
};

//bool operator==(side sa, side sb)
//{
//	if ((sa.a == sb.a) && (sa.b == sb.b))
//		return true;
//	else
//		return false;
//}

bool operator<(side sa, side sb)
{
	if (sa.a != sb.a)
		return (sa.a < sb.a);
	else
		return (sa.b < sb.b);
}



double getGeod(surface surface1, surface surface2)	//测地线距离
{
	double distx = surface1.centerx - surface2.centerx;
	double disty = surface1.centery - surface2.centery;
	double distz = surface1.centerz - surface2.centerz;

	return sqrt(distx * distx + disty * disty + distz * distz);
}

double getAngd(surface surface1, surface surface2)	//角距离，若是凹面，取eta=1，凸面取较小的0.1
{
	//先算两个法向量的模
	double mod1 = sqrt(surface1.nmvx * surface1.nmvx + surface1.nmvy * surface1.nmvy + surface1.nmvz * surface1.nmvz);
	double mod2 = sqrt(surface2.nmvx * surface2.nmvx + surface2.nmvy * surface2.nmvy + surface2.nmvz * surface2.nmvz);
	//两个法向量的内积
	double mul = surface1.nmvx * surface2.nmvx + surface1.nmvy * surface2.nmvy + surface1.nmvz * surface2.nmvz;
	double cos = mul / (mod1 * mod2);	//两个法向量的夹角的余弦

	//然后判断两个面的凹凸
	double ctcx = surface2.centerx - surface1.centerx;
	double ctcy = surface2.centery - surface1.centery;
	double ctcz = surface2.centerz - surface1.centerz;	//向量ctc是由面片1的质心指向面片2的质心
	//如果两个面是凸面，s1的法向量与ctc成钝角,内积为负，反之则成锐角,内积为正

	mul = ctcx * surface1.nmvx + ctcy * surface1.nmvy + ctcz * surface1.nmvz;

	double ratio = 1.0;
	if (mul < 0)
		ratio = eta;

	return ratio * (1 - cos);

}

double getDist(surface surface1, surface surface2) //总距离(这个总距离没有除以平均距离)
{
	return theta * getGeod(surface1, surface2) + (1 - theta) * getAngd(surface1, surface2);
}
