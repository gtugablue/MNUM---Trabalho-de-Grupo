#ifndef MONOCOMPARTIMENTAL_H
#define MONOCOMPARTIMENTAL_H

class Monocompartimental
{
private:
	double (*D)(double t);
public:
	Monocompartimental(double (*D)(double t));
};

#endif