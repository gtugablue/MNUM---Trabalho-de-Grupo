#ifndef _DOSAGEM_H
#define _DOSAGEM_H

#include <vector>

struct Toma
{
	double quantidade;
	double inicio;
};

class Dosagem
{
private:
	std::vector<Toma> tomas;
public:
	Dosagem();
	double operator()(double t) const;
};

#endif