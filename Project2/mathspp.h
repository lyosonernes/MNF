#pragma once
#include<vector>
class mathspp
{
public:
	mathspp();
	~mathspp();
	std::vector < std::vector<double>> DecCholesky(std::vector < std::vector<double>> A);
	double myErfInv2(double x);
};

