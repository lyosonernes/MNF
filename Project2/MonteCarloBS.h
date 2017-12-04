#pragma once
#include <vector>
#include "BlackScholes.h"

class MonteCarloBS
{
public:
	MonteCarloBS(double K, double r, double sig, double T,double Sz);
	std::pair <double,double> CalcCall(int N, int method, std::vector<double> Ti, double M, int VR, double c, int VRC);
	~MonteCarloBS();

private:
	double calcPayoffEU();
	std::pair<double, double> calcPayoffEUa();
	double calcPayoffAS(std::vector<double> St);
	double calcPayoffLB(std::vector<double> St);
	double calcPayoffUO(std::vector<double> St, double M);
	double calcPayoffLBSimCond(double ST);
	double calcPayoffUOSimCond(double ST, double M);
	double calcX(int method, std::vector<double> Ti,double M);
	std::pair<double, double> calcXa(int method, std::vector<double> Ti, double M);
	double maxCondBrow(double ST);
	double mbs_K;
	double mbs_r;
	double mbs_sig;
	double mbs_T;
	double mbs_Sz;
	BlackScholes mbs_bs;
	
};

