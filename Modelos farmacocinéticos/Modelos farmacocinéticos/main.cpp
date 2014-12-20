#include <iostream>
#include "Zeros.h"
#include <iomanip>
#include <time.h>
#include <fstream>
#include "Dosagem.h"

#define HORAS_POR_DIA			24
#define MINUTOS_POR_HORA		60
#define DOSE_DIARIA				500.0
#define NUM_HORAS_2A_TOMA		5.0 // Entre 4 e 6
#define VOL_APAR_PLASMA			3920.0
#define CONST_CIN_ELIM_TOTAL	0.001925
#define T_MAX					72.0

#define ERRO					1.0e-20

using namespace std;

double fKa(double Ka)
{
	return Ka * exp(-Ka * T_MAX) - CONST_CIN_ELIM_TOTAL * exp(-CONST_CIN_ELIM_TOTAL * T_MAX);
}

double gKa(double Ka)
{
	return (CONST_CIN_ELIM_TOTAL * T_MAX + log(Ka / CONST_CIN_ELIM_TOTAL)) / T_MAX;
}

double gKa2(double Ka)
{
	return CONST_CIN_ELIM_TOTAL * exp(Ka * T_MAX - CONST_CIN_ELIM_TOTAL * T_MAX);
}

double dfKa(double Ka)
{
	return exp(-Ka * T_MAX) - Ka * T_MAX * exp(-Ka * T_MAX);
}

int main()
{
	clock_t tStart;
	cout << setprecision(30);

	Dosagem D;
	ofstream outfile("dose.txt");
	for (double t = 0; t < 13 * 24 * 60; t += 0.5)
	{
		outfile << t << "\t" << D(t) << endl;
	}

	tStart = clock();
	cout << "---------- Metodo da bissecao ---------" << endl;
	cout << "Ka #1: " << metodoBissecao(fKa, 0, 0.01, ERRO) << endl;
	cout << "Ka #2: " << metodoBissecao(fKa, 0.01, 0.1, ERRO) << endl;
	cout << "Tempo de calculo: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;

	cout << endl;

	tStart = clock();
	cout << "---------- Metodo da corda ---------" << endl;
	cout << "Ka #1: " << metodoCorda(fKa, 0, 0.01, ERRO) << endl;
	cout << "Ka #2: " << metodoCorda(fKa, 0.01, 0.1, ERRO) << endl;
	cout << "Tempo de calculo: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;

	cout << endl;

	tStart = clock();
	cout << "---------- Metodo de Newton ---------" << endl;
	cout << "Ka #1: " << metodoNewton(fKa, dfKa, 0, ERRO) << endl;
	cout << "Ka #2: " << metodoNewton(fKa, dfKa, 0.05, ERRO) << endl;
	cout << "Tempo de calculo: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;

	cout << endl;

	tStart = clock();
	cout << "---------- Metodo de Picard-Peano ---------" << endl;
	cout << "\t---------- Funcao G1 ---------" << endl;
	cout << "Ka #1: " << metodoPicardPeano(gKa, 0.8, ERRO) << endl;
	cout << "Ka #2: " << metodoPicardPeano(gKa, 0.001, ERRO) << endl;
	cout << endl;
	cout << "\t---------- Funcao G2 ---------" << endl;
	cout << "Ka #1: " << metodoPicardPeano(gKa2, 0.002, ERRO) << endl;
	cout << "Ka #2: " << metodoPicardPeano(gKa2, 0.001, ERRO) << endl;
	cout << "Tempo de calculo: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;
	
	return 0;
}