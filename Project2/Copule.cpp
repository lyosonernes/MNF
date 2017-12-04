#include "Copule.h"
#include "mathspp.h"
#include "time.h"
#include <math.h>
#include <random>



Copule::Copule()
{
}


Copule::~Copule()
{
}

std::vector<double> Copule::gaussien(double rho, int n)
{
	mathspp math = mathspp();
	std::vector <double> U(n);

	// Construction de la matrice A
	std::vector<std::vector<double>> A;
	A.resize(n, std::vector<double>(n, 0));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (j == i) {
				A[i][j] = 1;
			}
			else {
				A[i][j] = rho;
			}
		}
	}

	// D�composition de Cholesky
	std::vector<std::vector<double>> L;
	L = math.DecCholesky(A);

	// G�n�ration uniforme ind�pendantes
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_real_distribution<double> distribution(0, 1);
	for (int i = 0; i < n; i++) {
		U[i] = distribution(generator);
	}

	//G�n�rateur du vecteur de gaussiennes ind�pendantes
	std::vector<double> Z(n);
	for (int i = 0; i < n; i++) {
		Z[i] = sqrt(2)*math.myErfInv2(2 * U[i] - 1);
		}

	//Vecteur de gaussiennes corr�l�es selon rho
	std::vector<double> Y(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Y[i] += L[j][i] * Z[i];
		}
	}

	//Vecteur de marginales uniformes li�es selon le copule gaussien de corr�lation rho
	std::vector<double> S(n);
	for (int i = 0; i < n; i++) {
		S[i] = (1 + erf(Y[i] / sqrt(2))) / 2;
	}

	return S;
}
