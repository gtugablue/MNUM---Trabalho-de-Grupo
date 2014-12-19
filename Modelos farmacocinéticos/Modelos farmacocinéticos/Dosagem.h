class Dosagem
{
private:
	double tempoToma;
	double quantidade;
public:
	Dosagem(const double tempoToma, const double quantidade);
	double operator()(double t) const;
};