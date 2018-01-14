#include <iostream>
#include "BlackScholes.h"
#include "Solver.h"
#include "MonteCarloBS.h"
#include "PriceA.h"
#include "mathspp.h"
#include "Copule.h"
#include <vector>
#include <random>
#include <iostream>
#include <fstream>

// FTD : T = 5
// Lambda i 0.05 
// N = 5,10,...

using namespace std;

int main()
{
	double K = 100; //100
	double r = 0.05; //0.05
	double sig = 0.2;
	double T = 1;
	double Si = 100;
	double t = 0.5;
	double L = log(1000);
	int N = 100000;
	int M = 1000;

	const int length = 20;



	BlackScholes mod1 = BlackScholes(K, r, sig, T);
	Solver m1 = Solver(L, T, sig, r, K);
	MonteCarloBS mcbs1 = MonteCarloBS(K, r, sig, T, Si);
	PriceA ps1 = PriceA(Si, r, sig);


	//// Initialisation d'un vecteur pour calculer BS
	//std::vector < std::vector<double>> Ubs;
	//Ubs.resize(M, std::vector<double>(N, 0));
	//Ubs = mod1.CalculBSU(N, M, L);

	//// Initialisation d'un vecteur pour calculer EE
	//std::vector < std::vector<double>> Uee;
	//Uee.resize(M, std::vector<double>(N, 0));
	//Uee = m1.CalculEulerE(N, M);
	

	//// Initialisation d'un vecteur pour calculer EI	
	//std::vector < std::vector<double>> Uei;
	//Uei.resize(M, std::vector<double>(N, 0));
	//Uei = m1.CalculEulerI(N, M);


	//// Export Ubs in Csv
	//ofstream myfile;
	//myfile.open("bsM1000N1(.csv");
	//for (int i = 0; i < N; i++) {
	//	myfile << exp(-L + i*(2 * L / N));
	//	myfile << ";";
	//	for (int j = 0; j < M; j++) {
	//		myfile << Ubs[j][i];
	//		myfile << ";";
	//	}
	//	myfile << "\n";
	//}

	//// Export Uee in Csv
	//ofstream myfile2;
	//myfile2.open("eulerexpliciteM1000N15.csv");
	//for (int i = 0; i < N; i++) {
	//	myfile2 << exp(-L + i*(2 * L / N));
	//	myfile2 << ";";
	//	for (int j = 0; j < M; j++) {
	//		myfile2 << Uee[j][i];
	//		myfile2 << ";";
	//	}
	//	myfile2 << "\n";
	//}

	// Création du vecteur Ti

	std::vector <double> Ti;
	Ti.resize(length);
	for (int i = 0; i < length; i++)
	{
		double x = 1 / double(length);
		Ti[i] = T*i*x;
	}

	//std::vector<double> St;
	//St.resize(Ti.size());
	//St = ps1.CalculSt(Ti);


	//// Export montecarlo in csv
	//double c0 = - mcbs1.calcC0(1000, Ti);
	//ofstream myfile;
	//myfile.open("mcbsOptionASRVC0.csv");
	//for (int i = 10; i < 15; i++) {
	//	std::pair<double, double> Ev = mcbs1.CalcCall(i, 1, Ti, 0, 2, c0, 0);
	//	myfile << i;
	//	myfile << ";";
	//	myfile << std::get<0>(Ev);
	//	myfile << ";";
	//	myfile << std::get<1>(Ev);
	//	myfile << "\n";
	//}

	//// Export montecarlo in csv
	////double c1 = -mcbs1.calcC1(1000, Ti);
	//double c1 = -0.46;
	//ofstream myfile2;
	//myfile2.open("mcbsOptionASRVC1.csv");
	//for (int i = 10; i < 15; i++) {
	//	std::pair<double, double> Ev = mcbs1.CalcCall(i, 1, Ti, 0, 2, c1, 1);
	//	myfile2 << i;
	//	myfile2 << ";";
	//	myfile2 << std::get<0>(Ev);
	//	myfile2 << ";";
	//	myfile2 << std::get<1>(Ev);
	//	myfile2 << "\n";
	//}

	//// Export montecarlo in csv
	//double c2 = -mcbs1.calcC2(1000, Ti);
	//ofstream myfile3;
	//myfile3.open("mcbsOptionASRVC2.csv");
	//for (int i = 10; i < 15; i++) {
	//	std::pair<double, double> Ev = mcbs1.CalcCall(i, 1, Ti, 0, 2, c2, 2);
	//	myfile3 << i;
	//	myfile3 << ";";
	//	myfile3 << std::get<0>(Ev);
	//	myfile3 << ";";
	//	myfile3 << std::get<1>(Ev);
	//	myfile3 << "\n";
	//}

	//// Export montecarlo in csv
	//ofstream myfile2;
	//myfile2.open("mcbsOptionAS.csv");
	//for (int i = 10; i < 9001; i++) {
	//	std::pair<double, double> Ev = mcbs1.CalcCall(i, 1, Ti, 0, 0, 0, 0);
	//	myfile2 << i;
	//	myfile2 << ";";
	//	myfile2 << std::get<0>(Ev);
	//	myfile2 << ";";
	//	myfile2 << std::get<1>(Ev);
	//	myfile2 << "\n";
	//}

	//// Export montecarlo in csv
	//ofstream myfile3;
	//myfile3.open("mcbsOptionASRVant.csv");
	//for (int i = 10; i < 9001; i++) {
	//	std::pair<double, double> Ev = mcbs1.CalcCall(i, 1, Ti, 0, 1, 0, 0);
	//	myfile3 << i;
	//	myfile3 << ";";
	//	myfile3 << std::get<0>(Ev);
	//	myfile3 << ";";
	//	myfile3 << std::get<1>(Ev);
	//	myfile3 << "\n";
	//}

	//// Export montecarlo in csv
	//ofstream myfile4;
	//myfile4.open("mcbsOptionDI.csv");
	//for (int i = 10; i < 9001; i++) {
	//	std::pair<double, double> Ev = mcbs1.CalcCall(i, 3, Ti, 50, 0, 0, 0);
	//	myfile4 << i;
	//	myfile4 << ";";
	//	myfile4 << std::get<0>(Ev);
	//	myfile4 << ";";
	//	myfile4 << std::get<1>(Ev);
	//	myfile4 << "\n";
	//}

	//// Export montecarlo in csv
	//ofstream myfile5;
	//myfile5.open("mcbsOptionDIant.csv");
	//for (int i = 10; i < 9001; i++) {
	//	std::pair<double, double> Ev = mcbs1.CalcCall(i, 3, Ti, 50, 1, 0, 0);
	//	myfile5 << i;
	//	myfile5 << ";";
	//	myfile5 << std::get<0>(Ev);
	//	myfile5 << ";";
	//	myfile5 << std::get<1>(Ev);
	//	myfile5 << "\n";
	//}

	//// Export St in csv
	//ofstream myfile;
	//myfile.open("prixactif.csv");
	//for (int i = 0; i < length; i++) {
	//	double test = St[i];
	//	myfile << St[i];
	//	myfile << "\n";
	//}

	//mathspp math = mathspp();

	//std::vector<std::vector<double>> A;
	//A.resize(3, std::vector<double>(3, 0));

	//A[0][0] = 1;
	//A[0][1] = 0.2;
	//A[0][2] = 0.2;
	//A[1][0] = 0.2;
	//A[1][1] = 1;
	//A[1][2] = 0.2;
	//A[2][0] = 0.2;
	//A[2][1] = 0.2;
	//A[2][2] = 1;

	//std::vector<std::vector<double>> C;
	////C.resize(2, std::vector<double>(2, 0));

	//C = math.DecCholesky(A);

	//// Export copule gaussien in csv
	//Copule cop = Copule();
	//ofstream myfile;
	//myfile.open("copulegauss0.8n.csv");
	//for (int i = 0; i < 500; i++) {
	//	std::vector<double> Ev = cop.gaussien(-0.8, 2);
	//	myfile << Ev[0];
	//	myfile << ";";
	//	myfile << Ev[1];
	//	myfile << "\n";
	//}


	//double a = mod1.CalculBS(0, Si);
	//cout << "10 , 20, 50 , 100 , 200 , 500 , 1000, 2000" << endl;
	//cout << "error EM" << endl;
	//cout << m1.ErrorDiscXtEMf(10, Si, N) << endl;
	//cout << m1.ErrorDiscXtEMf(20, Si, N) << endl;
	//cout << m1.ErrorDiscXtEMf(50, Si, N) << endl;
	//cout << m1.ErrorDiscXtEMf(100, Si, N) << endl;
	//cout << m1.ErrorDiscXtEMf(200, Si, N) << endl;
	//cout << m1.ErrorDiscXtEMf(500, Si, N) << endl;
	//cout << m1.ErrorDiscXtEMf(1000, Si, N) << endl;
	//cout << m1.ErrorDiscXtEMf(2000, Si, N) << endl;
	//cout << "error Mils" << endl;
	//cout << m1.ErrorDiscXtMilsf(10, Si, N) << endl;
	//cout << m1.ErrorDiscXtMilsf(20, Si, N) << endl;
	//cout << m1.ErrorDiscXtMilsf(50, Si, N) << endl;
	//cout << m1.ErrorDiscXtMilsf(100, Si, N) << endl;
	//cout << m1.ErrorDiscXtMilsf(200, Si, N) << endl;
	//cout << m1.ErrorDiscXtMilsf(500, Si, N) << endl;
	//cout << m1.ErrorDiscXtMilsf(1000, Si, N) << endl;
	//cout << m1.ErrorDiscXtMilsf(2000, Si, N) << endl;
	double a = m1.MCHeston(10000, 500, Si, 0.04, 0.04, 1, 0.5, -0.70, K);
	cout << a << endl;
	cout << mod1.CalculBS(0, Si) << endl;
	cout << mod1.CalcVolImpli(Si, a, 0.05, K) << endl;
	cout << "attente" << endl;
	return 0;
}