class Monocompartimental
{
private:
	double (*D)(double t);
public:
	Monocompartimental(double (*D)(double t));
};