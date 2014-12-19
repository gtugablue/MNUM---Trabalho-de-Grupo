#ifndef _ZEROS_H
#define _ZEROS_H

double metodoBissecao(double f(double), double a, double b, double erro);

double metodoCorda(double f(double), double a, double b, double erro);

double metodoNewton(double f(double), double df(double), double x, double erro);

double metodoPicardPeano(double g(double), double x, double erro);

#endif