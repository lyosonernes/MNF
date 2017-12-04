#pragma once
#include <vector>
class Copule
{
public:
	Copule();
	~Copule();
	std::vector<double> gaussien(double rho, int n);
};

