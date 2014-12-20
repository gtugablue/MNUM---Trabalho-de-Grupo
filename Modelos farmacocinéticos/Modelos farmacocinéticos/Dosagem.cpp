#include "Dosagem.h"
#include <cmath>

using namespace std;

#define HORAS_POR_DIA		24
#define MINUTOS_POR_HORA	60
#define DOSE_PEQUENA		50
#define DOSE_GRANDE			100
#define TEMPO_TOMA_POR_MG	0.06

Dosagem::Dosagem()
{
	Toma toma;

	toma.inicio = 0;
	toma.quantidade = DOSE_GRANDE;
	tomas.push_back(toma);
	
	toma.inicio += 4 * MINUTOS_POR_HORA;
	toma.quantidade = DOSE_GRANDE;
	tomas.push_back(toma);
	
	toma.inicio += 4 * MINUTOS_POR_HORA;
	toma.quantidade = DOSE_PEQUENA;
	tomas.push_back(toma);

	toma.inicio += 4 * MINUTOS_POR_HORA;
	toma.quantidade = DOSE_GRANDE;
	tomas.push_back(toma);

	toma.inicio += 4 * MINUTOS_POR_HORA;
	toma.quantidade = DOSE_PEQUENA;
	tomas.push_back(toma);

	toma.inicio += 4 * MINUTOS_POR_HORA;
	toma.quantidade = DOSE_GRANDE;
	tomas.push_back(toma);
}

double Dosagem::operator()(double t) const
{
	if (t > 12 * HORAS_POR_DIA * MINUTOS_POR_HORA)
	{
		return 0;
	}

	t = fmod(t, HORAS_POR_DIA * MINUTOS_POR_HORA);

	size_t i;
	for (i = 0; i < tomas.size(); ++i)
	{
		if (tomas[i].inicio < t && t < tomas[i].inicio + tomas[i].quantidade * TEMPO_TOMA_POR_MG)
		{
			return 1 / TEMPO_TOMA_POR_MG;
		}
	}

	return 0;
}