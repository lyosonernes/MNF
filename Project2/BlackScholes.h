#pragma once
#include <vector>
#include <array>
#include <iterator>

class BlackScholes
{
	
public:
	BlackScholes();
	BlackScholes(double K, double r, double sig, double T);
	~BlackScholes();
	float CalculBS(double t, double Si);
	float CalculBS(double t, double Si, double sig);
	std::vector < std::vector<double>> CalculBSU(int N, int M, double L);
	double CalcVolImpli(double S0, double PC, double sig1, double sig2);

private:
	double Calculbs_d1(double K, double r, double sig, double T, double Si, double t);
	double Calculbs_d2(double K, double r, double sig, double T, double Si, double t);
	double bs_d1;
	double bs_d2;
	double bs_K;
	double bs_r;
	double bs_sig;
	double bs_T;
	double bs_Si;
};