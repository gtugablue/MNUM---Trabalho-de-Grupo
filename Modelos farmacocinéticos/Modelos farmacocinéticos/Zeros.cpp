#include "Zeros.h"
#include <cmath>

double metodoBissecao(double f(double), double a, double b, double erro)
{
	double xold;
	double x = (a + b) / 2;

	do
	{
		if (f(a) * f(x) > 0)
			a = x;
		else
			b = x;

		xold = x;
		x = (a + b) / 2;
	} while (abs(x - xold) > erro);

	return x;
}

double metodoCorda(double f(double), double a, double b, double erro)
{
	double xold;
	double x = (a * f(b) - b * f(a)) / (f(b) - f(a));

	do
	{
		if (f(a) * f(x) > 0)
			a = x;
		else
			b = x;

		xold = x;
		x = (a * f(b) - b * f(a)) / (f(b) - f(a));
	} while (abs(x - xold) > erro);

	return x;
}

double metodoNewton(double f(double), double df(double), double x, double erro)
{
	double xold;

	do
	{
		xold = x;
		x = xold - f(xold) / df(xold);
	} while (abs(x - xold) > erro);

	return x;
}

double metodoPicardPeano(double g(double), double x, double erro)
{
	double xold;

	do
	{
		xold = x;
		x = g(xold);
	} while (abs(x - xold) > erro);

	return x;
}