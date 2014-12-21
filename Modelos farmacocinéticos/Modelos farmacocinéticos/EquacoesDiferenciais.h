#ifndef EQUACOES_DIFERENCIAIS_H
#define EQUACOES_DIFERENCIAIS_H

#include <vector>

class Point2D
{
public:
	Point2D(double x, double y) : x(x), y(y) {};
	double x;
	double y;
};

std::vector<Point2D> metodoEuler(double f(double, double), double x, double xf, double y, int num_steps);
std::vector<Point2D> metodoRungaKutta2a(double f(double, double), double x, double xf, double y, int num_steps);
std::vector<Point2D> metodoRungaKutta4a(double f(double, double), double x, double xf, double y, int num_steps);


#endif