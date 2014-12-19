#include "Bicompartimental.h"

Bicompartimental::Bicompartimental(const double Ke, const double Vap, double tMax, const Dosagem &D):
Ket(Ke / Vap), Ka()
{

}

const double Bicompartimental::fKa(const double Ka, const double tMax, const double Ke)
{

}