#include <iostream>
#include "Zeros.h"
#include <iomanip>
#include <string>
#include <time.h>
#include <fstream>
#include "Dosagem.h"
#include "EquacoesDiferenciais.h"
#include "SistemaEquacoesDiferenciais.h"

#define HORAS_POR_DIA			24
#define MINUTOS_POR_HORA		60
#define VOL_APAR_PLASMA			3920.0
#define CONST_CIN_ELIM_TOTAL	(0.1155 / MINUTOS_POR_HORA)
#define CONST_CIN_ELIM			(CONST_CIN_ELIM_TOTAL * VOL_APAR_PLASMA)
#define T_MAX					(1.2 * MINUTOS_POR_HORA)
#define FIM_DOSAGEM				(13 * HORAS_POR_DIA * MINUTOS_POR_HORA)
#define CONST_ABSORCAO			0.04600627564637752
//#define CONST_ABSORCAO			2.76038

#define ERRO					1.0e-20
#define NUM_ITER				(FIM_DOSAGEM * 128)

using namespace std;

Dosagem D;

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
	return D(t)/VOL_APAR_PLASMA - CONST_CIN_ELIM_TOTAL * Cp;
}

double bc1(double t, double mi, double mp)
{
	return D(t) - CONST_ABSORCAO * mi;
}

double bc2(double t, double mi, double mp)
{
	return CONST_ABSORCAO * mi - CONST_CIN_ELIM_TOTAL * mp;
}

double calcularQC(double S[])
{
	return (S[1] - S[0]) / (S[2] - S[1]);
}

void modeloMonocompartimental(Dosagem D)
{
	cout << "------------------------------------" << endl;
	cout << "---- Modelo monocompartimental -----" << endl;
	cout << "------------------------------------" << endl << endl;

	clock_t tStart;
	double tempoCalculo;
	vector<Point2D> pontos;
	double S[3];

	// EULER
	tStart = clock();
	cout << "--------- Metodo de Euler ---------" << endl << endl;
	pontos = metodoEuler(mc, 0, FIM_DOSAGEM, 0, NUM_ITER);
	tempoCalculo = ((double)clock() - tStart) / CLOCKS_PER_SEC;
	for (size_t i = pontos.size() - 1; i < pontos.size(); ++i)
	{
		cout << "t: " << pontos[i].x << "\tCp: " << pontos[i].y << endl;
	}
	cout << "Tempo de calculo: " << tempoCalculo << " s" << endl;

	ofstream outfile1("monocompartimental_euler.txt");
	for (size_t i = 0; i < pontos.size(); ++i)
	{
		outfile1 << pontos[i].x << "\t" << pontos[i].y << endl;
	}
	outfile1.close();
	cout << endl;

	for (size_t i = 0; i < 3; ++i)
	{
		pontos = metodoEuler(mc, 0, FIM_DOSAGEM, 0, pow(2, i) * NUM_ITER);
		S[i] = pontos[pontos.size() - 1].y;
		cout << "S[" << i << "] = " << S[i] << endl;
	}

	cout << "QC: " << calcularQC(S) << endl;

	cout << endl;

	// RK2
	tStart = clock();
	cout << "--------- Metodo de Runge-Kutta 2a ordem ---------" << endl << endl;
	pontos = metodoRungaKutta2a(mc, 0, FIM_DOSAGEM, 0, NUM_ITER);
	tempoCalculo = ((double)clock() - tStart) / CLOCKS_PER_SEC;
	for (size_t i = pontos.size() - 1; i < pontos.size(); ++i)
	{
		cout << "t: " << pontos[i].x << "\tCp: " << pontos[i].y << endl;
	}
	cout << "Tempo de calculo: " << tempoCalculo << " s" << endl;

	ofstream outfile2("monocompartimental_rk2.txt");
	for (size_t i = 0; i < pontos.size(); ++i)
	{
		outfile2 << pontos[i].x << "\t" << pontos[i].y << endl;
	}
	outfile2.close();
	cout << endl;

	for (size_t i = 0; i < 3; ++i)
	{
		pontos = metodoRungaKutta2a(mc, 0, FIM_DOSAGEM, 0, pow(2, i) * NUM_ITER);
		S[i] = pontos[pontos.size() - 1].y;
		cout << "S[" << i << "] = " << S[i] << endl;
	}

	cout << "QC: " << calcularQC(S) << endl;

	cout << endl;

	// RK4
	tStart = clock();
	cout << "--------- Metodo de Runge-Kutta 4a ordem ---------" << endl << endl;
	pontos = metodoRungaKutta4a(mc, 0, FIM_DOSAGEM, 0, NUM_ITER);
	tempoCalculo = ((double)clock() - tStart) / CLOCKS_PER_SEC;
	for (size_t i = pontos.size() - 1; i < pontos.size(); ++i)
	{
		cout << "t: " << pontos[i].x << "\tCp: " << pontos[i].y << endl;
	}
	cout << "Tempo de calculo: " << tempoCalculo << " s" << endl;

	ofstream outfile3("monocompartimental_rk4.txt");
	for (size_t i = 0; i < pontos.size(); ++i)
	{
		outfile3 << pontos[i].x << "\t" << pontos[i].y << endl;
	}
	outfile3.close();
	cout << endl;

	for (size_t i = 0; i < 3; ++i)
	{
		pontos = metodoRungaKutta4a(mc, 0, FIM_DOSAGEM, 0, pow(2, i) * NUM_ITER);
		S[i] = pontos[pontos.size() - 1].y;
		cout << "S[" << i << "] = " << S[i] << endl;
	}

	cout << "QC: " << calcularQC(S) << endl;

	cout << endl;
}

double bicompartimentalCalcularKa(Dosagem D)
{
	clock_t tStart;
	cout << setprecision(30);

	ofstream outfile("dose.txt");
	for (double t = 0; t < 13 * 24 * 60; t += 0.5)
	{
		outfile << t << "\t" << D(t) << endl;
	}

	cout << "------------------------------------" << endl;
	cout << "----- Modelo bicompartimental ------" << endl;
	cout << "------------------------------------" << endl;
	cout << "---------- Calculo do Ka -----------" << endl;
	cout << "------------------------------------" << endl << endl;

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

	return CONST_ABSORCAO;
}

void bicompartimentalCalcularMassas(Dosagem D)
{
	cout << "------------------------------------" << endl;
	cout << "----- Modelo bicompartimental ------" << endl;
	cout << "------------------------------------" << endl;
	cout << "-------- Calculo de massas ---------" << endl;
	cout << "------------------------------------" << endl << endl;

	clock_t tStart;
	double tempoCalculo;
	vector<Point3D> pontos;
	double Smi[3], Smp[3];

	// EULER
	tStart = clock();
	cout << "--------- Metodo de Euler ---------" << endl << endl;
	pontos = metodoEulerSistema(bc1, bc2, 0, FIM_DOSAGEM, 0, 0, NUM_ITER);
	tempoCalculo = ((double)clock() - tStart) / CLOCKS_PER_SEC;
	for (size_t i = pontos.size() - 1; i < pontos.size(); ++i)
	{
		cout << "t: " << pontos[i].x << "\tmi: " << pontos[i].y << "\tmp: " << pontos[i].z << endl;
	}
	cout << "Tempo de calculo: " << tempoCalculo << " s" << endl;

	ofstream outfile1("bicompartimental_euler.txt");
	for (size_t i = 0; i < pontos.size(); ++i)
	{
		outfile1 << pontos[i].x << "\t" << pontos[i].y << "\t" << pontos[i].z << endl;
	}
	outfile1.close();
	cout << endl;

	for (size_t i = 0; i < 3; ++i)
	{
		pontos = metodoEulerSistema(bc1, bc2, 0, FIM_DOSAGEM, 0, 0, pow(2, i) * NUM_ITER);
		Smi[i] = pontos[pontos.size() - 1].y;
		Smp[i] = pontos[pontos.size() - 1].z;
		cout << "Smi[" << i << "] = " << Smi[i] << endl;
		cout << "Smp[" << i << "] = " << Smp[i] << endl;
	}

	cout << "QCmi: " << calcularQC(Smi) << endl;
	cout << "QCmp: " << calcularQC(Smp) << endl;

	cout << endl;

	// RK2
	tStart = clock();
	cout << "--------- Metodo de Runge-Kutta 2a ordem ---------" << endl << endl;
	pontos = metodoRungaKutta2aSistema(bc1, bc2, 0, FIM_DOSAGEM, 0, 0, NUM_ITER);
	tempoCalculo = ((double)clock() - tStart) / CLOCKS_PER_SEC;
	for (size_t i = pontos.size() - 1; i < pontos.size(); ++i)
	{
		cout << "t: " << pontos[i].x << "\tmi: " << pontos[i].y << "\tmp: " << pontos[i].z << endl;
	}
	cout << "Tempo de calculo: " << tempoCalculo << " s" << endl;

	ofstream outfile2("bicompartimental_rk2.txt");
	for (size_t i = 0; i < pontos.size(); ++i)
	{
		outfile2 << pontos[i].x << "\t" << pontos[i].y << "\t" << pontos[i].z << endl;
	}
	outfile2.close();
	cout << endl;

	for (size_t i = 0; i < 3; ++i)
	{
		pontos = metodoRungaKutta2aSistema(bc1, bc2, 0, FIM_DOSAGEM, 0, 0, pow(2, i) * NUM_ITER);
		Smi[i] = pontos[pontos.size() - 1].y;
		Smp[i] = pontos[pontos.size() - 1].z;
		cout << "Smi[" << i << "] = " << Smi[i] << endl;
		cout << "Smp[" << i << "] = " << Smp[i] << endl;
	}

	cout << "QCmi: " << calcularQC(Smi) << endl;
	cout << "QCmp: " << calcularQC(Smp) << endl;

	cout << endl;

	// RK4
	tStart = clock();
	cout << "--------- Metodo de Runge-Kutta 4a ordem ---------" << endl << endl;
	pontos = metodoRungaKutta4aSistema(bc1, bc2, 0, FIM_DOSAGEM, 0, 0, NUM_ITER);
	tempoCalculo = ((double)clock() - tStart) / CLOCKS_PER_SEC;
	for (size_t i = pontos.size() - 1; i < pontos.size(); ++i)
	{
		cout << "t: " << pontos[i].x << "\tmi: " << pontos[i].y << "\tmp: " << pontos[i].z << endl;
	}
	cout << "Tempo de calculo: " << tempoCalculo << " s" << endl;

	ofstream outfile3("bicompartimental_rk4.txt");
	for (size_t i = 0; i < pontos.size(); ++i)
	{
		outfile3 << pontos[i].x << "\t" << pontos[i].y << "\t" << pontos[i].z << endl;
	}
	outfile3.close();

	for (size_t i = 0; i < 3; ++i)
	{
		pontos = metodoRungaKutta4aSistema(bc1, bc2, 0, FIM_DOSAGEM, 0, 0, pow(2, i) * NUM_ITER);
		Smi[i] = pontos[pontos.size() - 1].y;
		Smp[i] = pontos[pontos.size() - 1].z;
		cout << "Smi[" << i << "] = " << Smi[i] << endl;
		cout << "Smp[" << i << "] = " << Smp[i] << endl;
	}

	cout << "QCmi: " << calcularQC(Smi) << endl;
	cout << "QCmp: " << calcularQC(Smp) << endl;

	cout << endl;
}

void modeloBicompartimental(Dosagem D)
{
	while (1)
	{
		clearScreen();
		cout << "------------------------------------" << endl;
		cout << "----- Modelo bicompartimental ------" << endl;
		cout << "------------------------------------" << endl << endl;

		cout << "1. Calcular Ka" << endl;
		cout << "2. Calcular massas" << endl;
		cout << "3. Sair" << endl;

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
			bicompartimentalCalcularKa(D);
			pause();
			break;
		case 2:
			clearScreen();
			bicompartimentalCalcularMassas(D);
			pause();
			break;
		case 3:
			return;
		default:
			continue;
		}
	}
}

int main()
{
	cout << setprecision(20);
	Dosagem D;
	while (1)
	{
		clearScreen();
		cout << "------------------------------------" << endl;
		cout << "----- Modelos farmacocineticos -----" << endl;
		cout << "------------------------------------" << endl << endl;

		cout << "1. Modelo monocompartimental" << endl;
		cout << "2. Modelo bicompartimental" << endl;
		cout << "3. Sair" << endl;

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
			modeloMonocompartimental(D);
			pause();
			break;
		case 2:
			clearScreen();
			modeloBicompartimental(D);
			pause();
			break;
		case 3:
			return 0;
		default:
			continue;
		}
	}
}