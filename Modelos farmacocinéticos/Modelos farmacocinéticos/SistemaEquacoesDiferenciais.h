#ifndef SISTEMAS_EQUACOES_DIFERENCIAIS_H
#define SISTEMAS_EQUACOES_DIFERENCIAIS_H

#include <vector>

class Point3D
{
public:
	Point3D(double x, double y, double z) : x(x), y(y), z(z) {};
	double x;
	double y;
	double z;
};

std::vector<Point3D> metodoEulerSistema(double f1(double, double, double), double f2(double, double, double), double x, double xf, double y, double z, int num_steps);
std::vector<Point3D> metodoRungaKutta2aSistema(double f1(double, double, double), double f2(double, double, double), double x, double xf, double y, double z, int num_steps);
std::vector<Point3D> metodoRungaKutta4aSistema(double f1(double, double, double), double f2(double, double, double), double x, double xf, double y, double z, int num_steps);

#endif