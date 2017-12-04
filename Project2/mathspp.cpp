#include "mathspp.h"
#define _USE_MATH_DEFINES
#include <math.h>



mathspp::mathspp()
{
}


mathspp::~mathspp()
{
}

std::vector<std::vector<double>> mathspp::DecCholesky(std::vector<std::vector<double>> A)
{
	int n = A.size();
	std::vector<std::vector<double>> L;
	L.resize(n, std::vector<double>(n, 0));
	double a = 0;
	double b = 0;
	double c = 0;

	L[0][0] = sqrt(A[0][0]);
	for (int j = 1; j <= n - 1; j++) // Calcul de la première colonne
		L[j][0] = A[j][0] / L[0][0];
	for (int i = 1; i < n-1; i++)
	{
		for (int k = 0; k < i; k++)
			a += pow(L[i][k], 2);
		L[i][i] = sqrt(A[i][i] - a);
		for (int j =i + 1; j <n; j++)
		{
			for (int k = 0; k < i; k++)
				b += L[j][k] * L[i][k];
			L[j][i] = (A[j][i] - b) / L[i][i];
		}
	}
	for (int k = 0; k <n - 1; k++)
		c += pow(L[n - 1][k], 2);
	L[n - 1][n - 1] = sqrt(A[n - 1][n - 1] - c);
		return L;

}

double mathspp::myErfInv2(double x) { // CODE RECUPERE SUR INTERNET
	float tt1, tt2, lnx, sgn;
	sgn = (x < 0) ? -1.0 : 1.0;

	x = (1 - x)*(1 + x);        // x = 1 - x*x;
	lnx = logf(x);

	tt1 = 2 / (M_PI*0.147) + 0.5 * lnx;
	tt2 = 1 / (0.147) * lnx;

	return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}