#include "Solver.h"
#include <iterator>
#include <vector>
#include <math.h>
#include <algorithm> //std::max


Solver::Solver(double L, double T, double sig, double r,double K)
{
	S_Long = L;
	S_T = T;
	S_sig = sig;
	S_r = r;
	S_K = K;
}

Solver::~Solver()
{
}

std::vector < std::vector<double>> Solver::CalculEulerE(int N, int M)
{
	// Définition des pas
	double dx = (2 * S_Long) / N;
	double dt = S_T / M;

	// Création du vecteur U

	std::vector < std::vector<double>> U;
	U.resize(M, std::vector<double>(N, 0));

	// Création des arrays de discrétisation de l'espace et du temps

	double *y = new double[N];
	for (int i = 0; i <	 N; i++)
	{
		y[i] = -S_Long + i*dx;
	}

	double *t = new double[M];
	for (int j = 0; j < M; j++)
	{
		t[j] = j*dt;
	}


	// Initialisation du vecteur U (U0)

	for (int i = 0; i < N; i++)
	{
		U[0][i] = std::max(0.0, exp(y[i]) - S_K);
		//U[0][i] = 1;
	}

	// Définition des coefficients
	double ui = -(S_sig*S_sig / (dx*dx)) - S_r;
	double ui_p1 = (S_sig*S_sig / (2 * dx*dx)) + (S_r - S_sig*S_sig / 2)/(2*dx);
	double ui_m1 = (S_sig*S_sig / (2 * dx*dx)) - (S_r - S_sig*S_sig / 2)/(2*dx);

	for (int j = 0; j<M-1; j++)
	{
		if (j != 0) 
		{
			U[j][0] = 0; // cas du bord -L
			U[j][N - 1] = exp(y[N - 1]) - S_K*exp(-S_r*t[j]); //cas du bord L
		}


		for (int i = 1; i < N-1; i++)
		{
			U[j + 1][i] = U[j][i] + dt*(ui*U[j][i] + ui_p1*U[j][i + 1] + ui_m1*U[j][i - 1]);
		}

	}


	return U;
}

std::vector < std::vector<double>> Solver::CalculEulerI(int N, int M)
{
	// DÈfinition des pas
	double dx = 2 * S_Long / N;
	double dt = S_T / M;

	// CrÈation du vecteur U

	std::vector < std::vector<double>> U;
	U.resize(M, std::vector<double>(N, 0));

	// CrÈation des arays de discrÈtisation de l'espace et du temps

	double *y = new double[N];
	for (int i = 0; i <= N; i++)
	{
		y[i] = -S_Long + i*dx;
	};

	double *t = new double[M];
	for (int j = 0; j <= M; j++)
	{
		t[j] = j*dt;
	};


	// Initialisation du vecteur U (U0)

	for (int i = 0; i <= N; i++)
	{
		U[0][i] = std::max(0.0, exp(y[i]) - S_K);
	}

	// Definition des coefficients de A
	double ui = -(S_sig*S_sig / (dx*dx)) - S_r;
	double ui_p1 = (S_sig*S_sig / (2 * dx*dx)) + (S_r - S_sig*S_sig / 2) / (2 * dx);
	double ui_m1 = (S_sig*S_sig / (2 * dx*dx)) - (S_r - S_sig*S_sig / 2) / (2 * dx);

	//Definition des coefficients de ( I - DT A ) selon methode de Thomas
	double a = -dt*ui_m1;
	double b = 1 - dt*ui;
	double c = -dt*ui_p1;

	// Implementation algorithme de Thomas
	double *c_p = new double[N];
	double *d = new double[N];
	double *d_p = new double[N];

	for (int j = 0; j <= M - 2; j++)
	{
		// Calcul du vecteur d
		for (int i = 0; i <= N - 1; i++)
		{
			d[i] = U[j][i];
		}

		// Calcul des valeurs c'0 et d'0
		c_p[0] = c / b;
		d_p[0] = d[0] / b;


		// Calcul des valeur c'i
		for (int i = 1; i <= N - 2; i++)
		{
			c_p[i] = c / (b - a*c_p[i - 1]);
		}


		// Calcul des valeurs d'i
		for (int i = 1; i <= N - 1; i++)
		{
			d_p[i] = (d[i] - a*d_p[i - 1]) / (b - a*c_p[i - 1]);
		}

		// Condiditions au bord en -L et en L
		U[j][0] = 0;
		U[j][N - 1] = exp(y[N - 1]) - S_K*exp(-S_r*t[j]);


		// Formule de Thomas pour calculer de U de manière descendante
		// en fonction des vecteurs c' et d'

		U[j + 1][N - 1] = d_p[N - 1];
		for (int i = N - 2; i >= 0; i--)
		{
			U[j + 1][i] = d_p[i] - c_p[i] * U[i + 1][j + 1];
		}
	}
	return U;
}