#include <iostream>
#include "Zeros.h"
#include <iomanip>
#include <string>
#include <time.h>
#include <fstream>
#include "Dosagem.h"
#include "Monocompartimental.h"
#include "EquacoesDiferenciais.h"

#define HORAS_POR_DIA			24
#define MINUTOS_POR_HORA		60
#define DOSE_DIARIA				500.0
#define NUM_HORAS_2A_TOMA		5.0 // Entre 4 e 6
#define VOL_APAR_PLASMA			3920.0
#define CONST_CIN_ELIM_TOTAL	0.001925
#define T_MAX					72.0
#define FIM_DOSAGEM				13 * 24 * 60

#define ERRO					1.0e-20

using namespace std;

void clearScreen()
{
	system("cls");
}

void pause()
{
	system("pause");
}

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

double mc(double t, double Cp)
{
	Dosagem D;
	return (D(t) - CONST_CIN_ELIM_TOTAL * Cp) / VOL_APAR_PLASMA;
}

void calcularKa(Dosagem D)
{
	clock_t tStart;
	cout << setprecision(30);

	ofstream outfile("dose.txt");
	for (double t = 0; t < 13 * 24 * 60; t += 0.5)
	{
		outfile << t << "\t" << D(t) << endl;
	}

	cout << "CALCULO DO Ka" << endl;

	tStart = clock();
	cout << "---------- Metodo da bissecao ---------" << endl;
	cout << "Ka #1: " << metodoBissecao(fKa, 0, 0.01, ERRO) << endl;
	cout << "Ka #2: " << metodoBissecao(fKa, 0.01, 0.1, ERRO) << endl;
	cout << "Tempo de calculo: " << ((double)clock() - tStart) / CLOCKS_PER_SEC << " s" << endl;

	cout << endl;

	tStart = clock();
	cout << "---------- Metodo da corda ---------" << endl;
	cout << "Ka #1: " << metodoCorda(fKa, 0, 0.01, ERRO) << endl;
	cout << "Ka #2: " << metodoCorda(fKa, 0.01, 0.1, ERRO) << endl;
	cout << "Tempo de calculo: " << ((double)clock() - tStart) / CLOCKS_PER_SEC << " s" << endl;

	cout << endl;

	tStart = clock();
	cout << "---------- Metodo de Newton ---------" << endl;
	cout << "Ka #1: " << metodoNewton(fKa, dfKa, 0, ERRO) << endl;
	cout << "Ka #2: " << metodoNewton(fKa, dfKa, 0.05, ERRO) << endl;
	cout << "Tempo de calculo: " << ((double)clock() - tStart) / CLOCKS_PER_SEC << " s" << endl;

	cout << endl;

	tStart = clock();
	cout << "---------- Metodo de Picard-Peano ---------" << endl;
	cout << "\t---------- Funcao G1 ---------" << endl;
	cout << "Ka #1: " << metodoPicardPeano(gKa, 0.001925, ERRO) << endl;
	cout << "Ka #2: " << metodoPicardPeano(gKa, 0.046, ERRO) << endl;
	cout << "\t---------- Funcao G2 ---------" << endl;
	cout << "Ka #1: " << metodoPicardPeano(gKa2, 0.001925, ERRO) << endl;
	cout << "Ka #2: " << metodoPicardPeano(gKa2, 0.046, ERRO) << endl;
	cout << "Tempo de calculo: " << ((double)clock() - tStart) / CLOCKS_PER_SEC << " s" << endl;
}

void modeloMonocompartimental(Dosagem D)
{
	cout << "MODELO MONOCOMPARTIMENTAL" << endl;
	Monocompartimental monoCompartimental(D, VOL_APAR_PLASMA, CONST_CIN_ELIM_TOTAL * VOL_APAR_PLASMA);

	clock_t tStart;
	
	tStart = clock();
	cout << "--------- Metodo de Euler ---------" << endl;
	vector<Point> pontos = metodoEuler(mc, 0, FIM_DOSAGEM, 0, 1000);
	double tempoCalculo = ((double)clock() - tStart) / CLOCKS_PER_SEC;
	/*for (size_t i = 0; i < pontos.size(); ++i)
	{
	cout << "t: " << pontos[i].x << "\tCp: " << pontos[i].y << endl;
	}*/
	cout << "Tempo de calculo: " << tempoCalculo << " s" << endl;

	ofstream outfile("monocompartimental_euler.txt");
	for (size_t i = 0; i < pontos.size(); ++i)
	{
		outfile << pontos[i].x << "\t" << pontos[i].y << endl;
	}
}

int main()
{
	Dosagem D;
	while (1)
	{
		clearScreen();
		cout << "------------------------------------" << endl;
		cout << "----- Modelos farmacocineticos -----" << endl;
		cout << "------------------------------------" << endl;

		cout << endl;

		cout << "1. Calcular Ka" << endl;
		cout << "2. Modelo monocompartimental" << endl;
		cout << "3. Modelo bicompartimental" << endl;
		cout << "4. Sair" << endl;

		cout << endl;

		cout << "Insira uma opcao: " << endl;

		string option;
		getline(cin, option);
		unsigned optNumber;
		try
		{
			optNumber = stoi(option);
		}
		catch (...)
		{
			continue;
		}
		switch (optNumber)
		{
		case 1:
			clearScreen();
			calcularKa(D);
			pause();
			break;
		case 2:
			clearScreen();
			modeloMonocompartimental(D);
			pause();
			break;
		case 3:
		case 4:
			return 0;
		default:
			continue;
		}
	}
}