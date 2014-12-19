#include "Dosagem.h"
#include <cmath>

#define HORAS_POR_DIA		24
#define MINUTOS_POR_HORA	60

Dosagem::Dosagem(const double tempoToma, const double quantidade) :
tempoToma(tempoToma), quantidade(quantidade)
{

}

double Dosagem::operator()(double t) const
{
	t = fmod(t, HORAS_POR_DIA);

	if (t < MINUTOS_POR_HORA)
	{
		return 100 / MINUTOS_POR_HORA;
	}
	else if (t < tempoToma * 60)
	{
		return (quantidade - 100) / ((tempoToma - 1) * MINUTOS_POR_HORA);
	}

	return 0;
}