#include <iostream>
#include "Monocompartimental.h"

#define HORAS_POR_DIA			24
#define MINUTOS_POR_HORA		60
#define DOSE_DIARIA				500.0
#define NUM_HORAS_2A_TOMA		5.0 // Entre 4 e 6
#define VOL_APAR_PLASMA			3920.0
#define CONST_CIN_ELIM_TOTAL	0.1155

using namespace std;

double D(double t)
{
	t = fmod(t, HORAS_POR_DIA);

	if (t < MINUTOS_POR_HORA)
	{
		return 100 / MINUTOS_POR_HORA;
	}
	else if (t < (NUM_HORAS_2A_TOMA + 1) * 60)
	{
		return (DOSE_DIARIA - 100) / (NUM_HORAS_2A_TOMA * MINUTOS_POR_HORA);
	}

	return 0;
}

int main()
{
	return 0;
}