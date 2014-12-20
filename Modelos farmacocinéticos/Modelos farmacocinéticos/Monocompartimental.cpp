#include "Monocompartimental.h"

Monocompartimental::Monocompartimental(Dosagem D, double Vap, double Ke) : D(D), Vap(Vap), Ke(Ke)
{

}

double Monocompartimental::operator()(double t, double Cp)
{
	return (D(t) - Ke * Cp) / Vap;
}