#ifndef EQUACOES_DIFERENCIAIS_H
#define EQUACOES_DIFERENCIAIS_H

double metodoEuler(double f(double, double), double x, double xf, double y, int num_steps);
double metodoRungaKutta2a(double f(double, double), double x, double xf, double y, int num_steps);
double metodoRungaKutta4a(double f(double, double), double x, double xf, double y, int num_steps);

#endif