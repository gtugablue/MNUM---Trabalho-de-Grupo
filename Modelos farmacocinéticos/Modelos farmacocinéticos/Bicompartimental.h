#ifndef BICOMPARTIMENTAL_H
#define BICOMPARTIMENTAL_H

#include "Dosagem.h"

class Bicompartimental
{
private:
	static const double fKa(const double tMax, const double Ke);
	const double Ket;
	const double Ka;
	const Dosagem D;
public:
	Bicompartimental(const double Ket, const double Vap, double tMax, const Dosagem &D);
};

#endif