#include "PriceA.h"
#include <vector>
#include <array>
#include <random>

PriceA::PriceA(double S0, double r, double sig)
{
	ps_sz = S0;
	ps_r = r;
	ps_sig= sig;
}

PriceA::~PriceA()
{
}

std::vector<double> PriceA::CalculSt(std::vector<double> Ti)
{
	int length = Ti.size();
	std::vector<double> St;
	St.resize(length);

	std::random_device rd;
	std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0, 1);

	St[0] = ps_sz;

	for (int i = 1; i < length; i++) {
		double Z = distribution(generator);
		St[i] = St[i - 1] * exp((ps_r - ps_sig*ps_sig / 2)*(Ti[i] - Ti[i - 1]) + ps_sig*sqrt(Ti[i] - Ti[i - 1])*Z);
	}

	return St;
}

std::pair<std::vector<double>, std::vector<double>> PriceA::CalculSta(std::vector<double> Ti)
{
	std::pair<std::vector<double>, std::vector<double>> Sta;
	int length = Ti.size();
	std::get<0>(Sta).resize(length);
	std::get<1>(Sta).resize(length);

	std::random_device rd;
	std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0, 1);

	std::get<0>(Sta)[0] = ps_sz;
	std::get<1>(Sta)[0] = ps_sz;

	for (int i = 1; i < length; i++) {
		double Z = distribution(generator);
		std::get<0>(Sta)[i] = std::get<0>(Sta)[i - 1] * exp((ps_r - ps_sig*ps_sig / 2)*(Ti[i] - Ti[i - 1]) + ps_sig*sqrt(Ti[i] - Ti[i - 1])*Z);
		std::get<1>(Sta)[i] = std::get<1>(Sta)[i - 1] * exp((ps_r - ps_sig*ps_sig / 2)*(Ti[i] - Ti[i - 1]) - ps_sig*sqrt(Ti[i] - Ti[i - 1])*Z);
	}

	return Sta;
}
