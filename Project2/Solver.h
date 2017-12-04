#pragma once
#include <vector>
#include <array>
#include <iterator>

class Solver
{
public:
	Solver(double L, double T, double sig, double r, double K);
	~Solver();
	std::vector < std::vector<double>> CalculEulerE(int N, int M);
	std::vector < std::vector<double>> CalculEulerI(int N, int M);

private:
	double S_Long;
	double S_T;
	double S_sig;
	double S_r;
	double S_K;

};

