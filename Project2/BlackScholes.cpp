#include "BlackScholes.h"
#define _USE_MATH_DEFINES
#include <math.h>




BlackScholes::BlackScholes()
{
}

BlackScholes::BlackScholes(double K, double r, double sig, double T)
{
	//Initialisation du mod�le avec les param�tres K,r,sig,T
	bs_K = K;
	bs_r = r;
	bs_sig = sig;
	bs_T = T;
}

BlackScholes::~BlackScholes()
{
}


float BlackScholes::CalculBS(double t, double Si)
{

	// si le prix de l'actif est nul, la valeur de l'option est nulle
	if (Si < 1.e-14) return 0;

	// Si sigma est null, nous devons selon le signe de d1 et d2 retourner 0
	if (bs_sig < 1.e-14)
	{
		if (Si < bs_K*exp(-bs_r*(bs_T - t)))return 0.;
		else return Si - bs_K*exp(-bs_r*(bs_T - t));
	}

	// Si on demande le prix au temps du strike, selon le signe de ST-K, on retourne le strike ou 0.
	if (fabs(bs_T - t)<1.e-14)
	{
		if (Si < bs_K)return 0.;
		else return Si - bs_K;
	}

	//Calcul de d1 et d2
	bs_d1 = Calculbs_d1(bs_K, bs_r, bs_sig, bs_T, Si,t);
	bs_d2 = Calculbs_d2(bs_K, bs_r, bs_sig, bs_T, Si,t);
	double a = (1 + erf(bs_d1 / sqrt(2))) / 2;
	double b = (1 + erf(bs_d2 / sqrt(2))) / 2;
	return Si*a - bs_K*exp(-bs_r*(bs_T-t))*b;


}

float BlackScholes::CalculBS(double t, double Si,double sig, double K)
{

	// si le prix de l'actif est nul, la valeur de l'option est nulle
	if (Si < 1.e-14) return 0;

	// Si sigma est null, nous devons selon le signe de d1 et d2 retourner 0
	if (sig < 1.e-14)
	{
		if (Si < K*exp(-bs_r*(bs_T - t)))return 0.;
		else return Si - K*exp(-bs_r*(bs_T - t));
	}

	// Si on demande le prix au temps du strike, selon le signe de ST-K, on retourne le strike ou 0.
	if (fabs(bs_T - t)<1.e-14)
	{
		if (Si < K)return 0.;
		else return Si - K;
	}

	//Calcul de d1 et d2
	bs_d1 = Calculbs_d1(K, bs_r, sig, bs_T, Si, t);
	bs_d2 = Calculbs_d2(K, bs_r, sig, bs_T, Si, t);
	double a = (1 + erf(bs_d1 / sqrt(2))) / 2;
	double b = (1 + erf(bs_d2 / sqrt(2))) / 2;
	return Si*a - K*exp(-bs_r*(bs_T - t))*b;


}

std::vector<std::vector<double>> BlackScholes::CalculBSU(int N, int M, double L)
{

	// D�finition des pas
	double dx = (2 * L) / N;
	double dt = bs_T / M;

	// Cr�ation du vecteur U

	std::vector < std::vector<double>> U;
	U.resize(M, std::vector<double>(N, 0));

	double *y = new double[N];
	for (int i = 0; i < N; i++)
	{
		y[i] = -L + i*dx;
	}

	for (int j = 0; j<M ; j++)
	{
		for (int i = 0; i < N ; i++)
		{
			U[j][i] = CalculBS(bs_T - j*dt, exp(y[i]));
		}

	}

	return U;
}

double BlackScholes::CalcVolImpli(double S0, double PC, double sig1, double K)
{
	double test= CalculBS(0, S0, sig1, K);
	double signe = 1;
	if (PC < signe) signe = -1 ;
	double xnp1 = sig1;
	double xn = 0;
	while (abs(xnp1 - xn) > 1.e-3) {
		xn = xnp1;
		double d1 = Calculbs_d1(K, bs_r, xn, bs_T, S0, 0);
		double vega = (exp(-pow(d1,2)/2)*S0*sqrt(bs_T))/ (sqrt(2 * M_PI));
		xnp1 = xn - signe*((CalculBS(0, S0, xn, K) - PC) / vega);
	}
	return xnp1;
}


//double BlackScholes::CalcVolImpli(double S0, double PC, double sig1, double sig2, double K)
//{
//	double a = sig1;
//	double b = sig2;
//	double c = (a + b) / 2;
//	double cbs = CalculBS(bs_T, S0, c, K);
//	while (abs(cbs - PC) > 1.e-3) {
//		if (PC > cbs) {
//			a = c;
//		}
//		else {
//			b = c;
//		}
//		c = (a + b) / 2;
//		cbs = CalculBS(0, S0, c,K);
//	}
//	return c;
//}

double BlackScholes::Calculbs_d1(double K, double r, double sig, double T,double Si,double t)
{
	return (log(Si / K) + (r + (sig*sig) / 2)*(T - t)) / (sig*sqrt(T - t));
}

double BlackScholes::Calculbs_d2(double K, double r, double sig, double T,double Si,double t)
{
	return (log(Si / K) + (r + (sig*sig) / 2)*(T - t)) / (sig*sqrt(T - t)) - sig*sqrt(T-t);
}






