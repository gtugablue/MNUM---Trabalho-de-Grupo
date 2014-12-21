#include "Zeros.h"
#include <cmath>
#include <iostream>

using namespace std;

double metodoBissecao(double f(double), double a, double b, double erro)
{
	double xold;
	double x = (a + b) / 2;
	size_t i = 0;
	do
	{
		if (f(a) * f(x) > 0)
			a = x;
		else
			b = x;

		xold = x;
		x = (a + b) / 2;
		++i;
		if (i >= 100) break;
	} while (abs(x - xold) > erro);

	cout << "Numero de iteracoes: " << i << endl;

	return x;
}

double metodoCorda(double f(double), double a, double b, double erro)
{
	double xold;
	double x = (a * f(b) - b * f(a)) / (f(b) - f(a));
	size_t i = 0;
	do
	{
		if (f(a) * f(x) > 0)
			a = x;
		else
			b = x;

		xold = x;
		x = (a * f(b) - b * f(a)) / (f(b) - f(a));
		++i;
		if (i >= 100) break;
	} while (abs(x - xold) > erro);

	cout << "Numero de iteracoes: " << i << endl;

	return x;
}

double metodoNewton(double f(double), double df(double), double x, double erro)
{
	double xold;
	size_t i = 0;
	do
	{
		xold = x;
		x = xold - f(xold) / df(xold);
		++i;
		if (i >= 100) break;
	} while (abs(x - xold) > erro);

	cout << "Numero de iteracoes: " << i << endl;

	return x;
}

double metodoPicardPeano(double g(double), double x, double erro)
{
	double xold;
	size_t i = 0;
	do
	{
		xold = x;
		x = g(xold);
		++i;
		if (i >= 100) break;
	} while (abs(x - xold) > erro);

	cout << "Numero de iteracoes: " << i << endl;

	return x;
}