#ifndef EQUACOES_DIFERENCIAIS_H
#define EQUACOES_DIFERENCIAIS_H
#include <vector>

class Point
{
public:
	Point(double x, double y) : x(x), y(y) {};
	double x;
	double y;
};

std::vector<Point> metodoEuler(double f(double, double), double x, double xf, double y, int num_steps);
std::vector<Point> metodoRungaKutta2a(double f(double, double), double x, double xf, double y, int num_steps);
std::vector<Point> metodoRungaKutta4a(double f(double, double), double x, double xf, double y, int num_steps);


#endif