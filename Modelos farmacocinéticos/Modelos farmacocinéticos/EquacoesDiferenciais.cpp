#include "EquacoesDiferenciais.h"
#include <cmath>

double metodoEuler(double f(double, double), double x, double xf, double y, int num_steps)
{
	double h = (xf - x) / num_steps;

	for (int i = 0; i < num_steps; ++i)
	{
		y += h * f(x, y);
		x += h;
	}

	return y;
}

double metodoRungaKutta2a(double f(double, double), double x, double xf, double y, int num_steps)
{
	double h = (xf - x) / num_steps;

	for (int i = 0; i < num_steps; ++i)
	{
		y += h * f(x + h / 2, y + h / 2 * f(x, y));
		x += h;
	}
	return y;
}

double metodoRungaKutta4a(double f(double, double), double x, double xf, double y, int num_steps)
{
	double h = (xf - x) / num_steps;
	double deltaY1, deltaY2, deltaY3, deltaY4;
	for (int i = 0; i < num_steps; ++i)
	{
		deltaY1 = h * f(x, y);
		deltaY2 = h * f(x + h / 2, y + deltaY1 / 2);
		deltaY3 = h * f(x + h / 2, y + deltaY2 / 2);
		deltaY4 = h * f(x + h, y + deltaY3);

		x += h;
		y += (1.0 / 6) * deltaY1 + (1.0 / 3) * deltaY2 + (1.0 / 3) * deltaY3 + (1.0 / 6) * deltaY4;
	}
	return y;
}