#include "SistemaEquacoesDiferenciais.h"

using namespace std;

vector<Point3D> metodoEulerSistema(double f1(double, double, double), double f2(double, double, double), double x, double xf, double y, double z, int num_steps)
{
	double h = (xf - x) / num_steps;
	vector<Point3D> res;
	res.push_back(Point3D(x, y, z));

	for (int i = 0; i < num_steps; ++i)
	{
		z += h * f2(x, y, z);
		y += h * f1(x, y, z);
		x += h;
		res.push_back(Point3D(x, y, z));
	}

	return res;
}

vector<Point3D> metodoRungaKutta2aSistema(double f1(double, double, double), double f2(double, double, double), double x, double xf, double y, double z, int num_steps)
{
	double h = (xf - x) / num_steps;
	double deltaY, deltaK;
	vector<Point3D> res;
	res.push_back(Point3D(x, y, z));

	for (int i = 0; i < num_steps; ++i)
	{
		deltaY = h * f1(x + h / 2, y + h / 2 * f1(x, y, z), z + h / 2 * f2(x, y, z));
		deltaK = h * f2(x + h / 2, y + h / 2 * f1(x, y, z), z + h / 2 * f2(x, y, z));
	
		x += h;
		y += deltaY;
		z += deltaK;
		res.push_back(Point3D(x, y, z));
	}
	return res;
}

vector<Point3D> metodoRungaKutta4aSistema(double f1(double, double, double), double f2(double, double, double), double x, double xf, double y, double z, int num_steps)
{
	double h = (xf - x) / num_steps;
	double deltaY1, deltaY2, deltaY3, deltaY4, deltaK1, deltaK2, deltaK3, deltaK4;
	vector<Point3D> res;
	res.push_back(Point3D(x, y, z));

	for (int i = 0; i < num_steps; ++i)
	{
		deltaY1 = h * f1(x, y, z);
		deltaK1 = h * f2(x, y, z);

		deltaY2 = h * f1(x + h / 2, y + deltaY1 / 2, z + deltaK1 / 2);
		deltaK2 = h * f2(x + h / 2, y + deltaY1 / 2, z + deltaK1 / 2);

		deltaY3 = h * f1(x + h / 2, y + deltaY2 / 2, z + deltaK2 / 2);
		deltaK3 = h * f2(x + h / 2, y + deltaY2 / 2, z + deltaK2 / 2);

		deltaY4 = h * f1(x + h, y + deltaY3, z + deltaK3);
		deltaK4 = h * f2(x + h, y + deltaY3, z + deltaK3);

		x += h;
		y += deltaY1 / 6 +  deltaY2 / 3 + deltaY3 / 3 + deltaY4 / 6;
		z += deltaK1 / 6 + deltaK2 / 3 + deltaK3 / 3 + deltaK4 / 6;

		res.push_back(Point3D(x, y, z));
	}
	return res;
}
