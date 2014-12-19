#include "EquacoesDiferenciais.h"
#include <cmath>

double metodoEuler(double(*f)(double x, double y), double xi, double xf, double y, double h)
{
	unsigned n = (xf - xi) / h;
	double newX = xi;
	for (size_t i = 0; i <= n; ++i)
	{
		y += h * f(newX, y);
		newX += h;
	}
	return pow(xi, 2) + pow(y, 2);
}

double metodoRungaKutta2a(double f(double x, double y), double xi, double xf, double y, double h)
{
	unsigned n = (xf - xi) / h;
	double newX = xi;
	for (size_t i = 0; i <= n; ++i)
	{
		y += h * f(newX + h / 2, y + h / 2 * f(newX, y));
		newX += h;
	}
	return y;
}

double metodoRungaKutta4a(double f(double x, double y), double xi, double xf, double y, double h)
{
	unsigned n = (xf - xi) / h;
	double deltaY1, deltaY2, deltaY3, deltaY4;
	for (size_t i = 0; i <= n; ++i)
	{
		deltaY1 = h * f(xi, y);
		deltaY2 = h * f(xi + h / 2, y + deltaY1 / 2);
		deltaY3 = h * f(xi + h / 2, y + deltaY2 / 2);
		deltaY4 = h * f(xi + h, y + deltaY3);
		xi += h;
		y += (1.0 / 6) * deltaY1 + (1.0 / 3) * deltaY2 + (1.0 / 3) * deltaY3 + (1.0 / 6) * deltaY4;
	}
	return y;
}