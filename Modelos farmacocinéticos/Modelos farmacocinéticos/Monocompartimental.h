#ifndef MONOCOMPARTIMENTAL_H
#define MONOCOMPARTIMENTAL_H

#include "Dosagem.h"

class Monocompartimental
{
private:
	Dosagem D;
	double Vap;
	double Ke;
public:
	Monocompartimental(Dosagem D, double Vap, double Ke);
	double operator()(double t, double Cp);
	//operator double(*)(double, double)() { return operator() }
};

/*double(*)(double, double) static_cast<double(*)(double, double)>(Monocompartimental monocompartimental)
{
	return monocompartimental.operator();
}*/

#endif