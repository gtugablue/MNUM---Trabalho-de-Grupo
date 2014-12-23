#include "EquacoesDiferenciais.h"
#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

vector<Point2D> metodoEuler(double f(double, double), double x, double xf, double y, int num_steps)
{
	double h = (xf - x) / num_steps;
	vector<Point2D> res;
	res.push_back(Point2D(x, y));
	for (int i = 0; i < num_steps; ++i)
	{
		y += h * f(x, y);
		x += h;
		res.push_back(Point2D(x, y));
	}

	return res;
}

vector<Point2D> metodoRungaKutta2a(double f(double, double), double x, double xf, double y, int num_steps)
{
	double h = (xf - x) / num_steps;
	vector<Point2D> res;
	res.push_back(Point2D(x, y));
	for (int i = 0; i < num_steps; ++i)
	{
		y += h * f(x + h / 2, y + h / 2 * f(x, y));
		x += h;
		res.push_back(Point2D(x, y));
	}
	return res;
}

vector<Point2D> metodoRungaKutta4a(double f(double, double), double x, double xf, double y, int num_steps)
{
	double h = (xf - x) / num_steps;
	double deltaY1, deltaY2, deltaY3, deltaY4;
	vector<Point2D> res;
	res.push_back(Point2D(x, y));
	for (int i = 0; i < num_steps; ++i)
	{
		deltaY1 = h * f(x, y);
		deltaY2 = h * f(x + h / 2, y + deltaY1 / 2);
		deltaY3 = h * f(x + h / 2, y + deltaY2 / 2);
		deltaY4 = h * f(x + h, y + deltaY3);

		y += (1.0 / 6) * deltaY1 + (1.0 / 3) * deltaY2 + (1.0 / 3) * deltaY3 + (1.0 / 6) * deltaY4;
		x += h;
		res.push_back(Point2D(x, y));
	}
	return res;
}