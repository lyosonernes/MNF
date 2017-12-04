#pragma once
#include <vector>
#include <array>
#include <iterator>

class PriceA
{
public:
	PriceA(double S0, double r, double sig);
	~PriceA();
	std::vector<double> CalculSt(std::vector<double> Ti);
	std::pair<std::vector<double>, std::vector<double>> CalculSta(std::vector<double> Ti);

private:
	double ps_sz;
	double ps_r;
	double ps_sig;

};

