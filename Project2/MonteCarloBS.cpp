#include "MonteCarloBS.h"
#include "BlackScholes.h"
#include "PriceA.h"
#include <random>
#include <algorithm>
#include <time.h>
#include <vector>
#include <numeric>



MonteCarloBS::MonteCarloBS(double K, double r, double sig, double T, double Sz)
{
	mbs_K = K;
	mbs_r = r;
	mbs_sig = sig;
	mbs_T = T;
	mbs_Sz = Sz;
	mbs_bs = BlackScholes(mbs_K, mbs_r, mbs_sig, mbs_T);
}

std::pair<double,double> MonteCarloBS::CalcCall(int N, int method = 0, std::vector<double> Ti = std::vector<double>() , double M = 0, int VR = 0, double c = 0, int VRC = 0)
{
	double Cm = 0;
	double Cmx = 0;
	double Cmy = 0;
	double Cmcarre = 0;
	double Cmxcarre = 0;
	double Cmycarre = 0;
	double sr = 0;
	double cetoile = 0;
	std::pair <double, double> Ev(0, 0);
	std::pair <double, double> Eva(0, 0);

	switch (VR) 
	{
		case 0: // Pas de réduction de variation
			for (int i = 0; i < N; i++)
			{
				sr = calcX(method, Ti, M);
				Cm = Cm + sr;
				Cmcarre = Cmcarre + sr*sr;
			}
			Cm = Cm / N;
			std::get<0>(Ev) = Cm;
			std::get<1>(Ev) = Cmcarre / N - Cm*Cm;
			return Ev ;
			break;
		case 1: // Réduction de variance : Antithétique
			for (int i = 0; i < N; i++)
			{
				Eva = calcXa(method, Ti, M);
				Cmx = Cmx + std::get<0>(Eva);
				Cmy = Cmy + std::get<1>(Eva);
				//Cm = Cm + std::get<0>(Eva)*std::get<1>(Eva);
				//Cmxcarre = Cmxcarre + std::get<0>(Eva)*std::get<0>(Eva);
				//Cmycarre = Cmycarre + std::get<1>(Eva)*std::get<1>(Eva);
				Cmcarre = Cmcarre + pow(std::get<0>(Eva),2) + pow(std::get<1>(Eva), 2);
			}
			Cmx = Cmx / N;
			Cmy = Cmy / N;
			//Cm = Cm / N;
			//Cmxcarre = Cmxcarre / N;
			//Cmycarre = Cmycarre / N;
			std::get<0>(Ev) = (Cmx + Cmy)/2 ;
			//std::get<1>(Ev) = ((Cmxcarre - Cmx*Cmx) + (Cmycarre - Cmy*Cmy) + 2 * (Cm - Cmx*Cmy)) / 4;
			std::get<1>(Ev) = (Cmcarre / (2*N)) - pow(std::get<0>(Ev), 2);
			return Ev;
			break;
		case 2 : // Réduction de variance : Contrôle
			for (int i = 0; i < N; i++)
			{
				switch (VRC)
				{
				case 0: // Variante Z = ST

					sr = calcX((method*(VR/2))*10 + VRC, Ti, c);
					break;
				case 1: // Variante Z = 1/N(somme Sti)
					sr = calcX((method*(VR/2))*10 + VRC, Ti, c);
					break;
				case 2: // Variante Z = exp(-rT)*(ST-K)+
					sr = calcX((method*(VR / 2))*10 + VRC, Ti, c);
					//sr = calcX(method, St, M) + c*(exp(-mbs_r*mbs_T)*std::max(0.0, 0.0) - mbsmod1.CalculBS(0, mbs_Sz));
					break;
				default:
					sr = 0.0;
				}
				Cm = Cm + sr;
				Cmcarre = Cmcarre + sr*sr;
			}
			Cm = Cm / N;
			std::get<0>(Ev) = Cm;
			std::get<1>(Ev) = Cmcarre / N - Cm*Cm; 
			return Ev;
			break;
		default :
			return Ev;
	}


}


MonteCarloBS::~MonteCarloBS()
{
}


double MonteCarloBS::calcPayoffEU()
{
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0, 1); // Peut être pas optimal niveau temps de calcul
	double Z = distribution(generator);
	return std::max(0.0,mbs_Sz*exp((mbs_r- mbs_sig*mbs_sig/2)*mbs_T+ mbs_sig*sqrt(mbs_T)*Z)-mbs_K);	
}

std::pair<double, double>  MonteCarloBS::calcPayoffEUa()
{
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0, 1); // Peut être pas optimal niveau temps de calcul
	std::pair <double, double> Ev(0, 0);
	double Z = distribution(generator);
	std::get<0>(Ev) = std::max(0.0, mbs_Sz*exp((mbs_r - mbs_sig*mbs_sig / 2)*mbs_T + mbs_sig*sqrt(mbs_T)*Z) - mbs_K);
	std::get<1>(Ev) = std::max(0.0, mbs_Sz*exp((mbs_r - mbs_sig*mbs_sig / 2)*mbs_T - mbs_sig*sqrt(mbs_T)*Z) - mbs_K);
	return Ev;
}


double MonteCarloBS::calcPayoffAS(std::vector<double> St)
{
	double stm = 0;
	int size = St.size();
	for (int i = 0; i < size; i++) {
		stm += St[i];
	}
	stm = stm / size;
	return std::max(0.0, stm - mbs_K);

}

double MonteCarloBS::calcPayoffLB(std::vector<double> St)
{
	double stm = 0;
	int size = sizeof(St);
	for (int i = 0; i < size; i++) {
		stm = std::max(stm, St[i]);
	}
	return std::max(0.0, stm - mbs_K);
}

double MonteCarloBS::calcPayoffUO(std::vector<double> St, double M)
{
	double stm = 0;
	int size = St.size();
	for (int i = 0; i < size; i++) {
		stm = std::min(stm, St[i]);
	}
	if (stm <= M) {
		stm = St[size-1];
	}
	return std::max(0.0, stm - mbs_K);
}

double MonteCarloBS::calcPayoffLBSimCond(double ST)
{
	double stmax = maxCondBrow(ST);
	return std::max(0.0, stmax - mbs_K);
}

double MonteCarloBS::calcPayoffUOSimCond(double ST, double M)
{
	double stmax = maxCondBrow(ST);
	if (stmax >= M) {
		stmax = 0;
	}
	return std::max(0.0, stmax - mbs_K);
}

double MonteCarloBS::calcX(int method, std::vector<double> Ti, double M = 0)
{
	PriceA ps1 = PriceA(mbs_Sz, mbs_r, mbs_sig);
	std::vector<double> St;
	double esp_z = 0;
	switch (method)
	{
		case 0 : //Call Européen
			return exp(-mbs_r*mbs_T)*calcPayoffEU();
			break;
		case 1 : //Call Asiatiaque
			St.resize(Ti.size());
			St = ps1.CalculSt(Ti);
			return exp(-mbs_r*mbs_T)*calcPayoffAS(St);
			break;
		case 10 : //Call asiatique VR condition 1 + VR*10*method + VRC
			St.resize(Ti.size());
			St = ps1.CalculSt(Ti);
			return exp(-mbs_r*mbs_T)*calcPayoffAS(St) + M*(St[Ti.size() -1] - mbs_Sz*exp(mbs_r*mbs_T));
			break;
		case 11 : //Call asiatique VR condition 2 + VR*10*method + VRC
			St.resize(Ti.size());
			St = ps1.CalculSt(Ti);
			for (int i = 0; i < Ti.size(); i++) {
				esp_z = esp_z + mbs_Sz*exp(mbs_r*Ti[i]);
			}
			esp_z = esp_z / Ti.size();
			return exp(-mbs_r*mbs_T)*calcPayoffAS(St) + M*(std::accumulate(St.begin(), St.end(), 0.0)/Ti.size() - esp_z);
			break;
		case 12 : //Call asiatique VR condition 3 + VR*10*method + VRC
			St.resize(Ti.size());
			St = ps1.CalculSt(Ti);
			return exp(-mbs_r*mbs_T)*calcPayoffAS(St) + M*(exp(-mbs_r*mbs_T)*std::max(0.0, St[Ti.size()-1] - mbs_K) - mbs_bs.CalculBS(0, mbs_Sz));
			break;
		case 2 : //Call Lookback discrète
			St.resize(Ti.size());
			St = ps1.CalculSt(Ti);
			return exp(-mbs_r*mbs_T)*calcPayoffLB(St);
			break;
		case 3 : //Call Up and out discrète
			St.resize(Ti.size());
			St = ps1.CalculSt(Ti);
			return exp(-mbs_r*mbs_T)*calcPayoffUO(St,M);
			break;
		default :
			return 0.0;
	}
}


std::pair<double, double> MonteCarloBS::calcXa(int method, std::vector<double> Ti,double M)
{
	PriceA ps1 = PriceA(mbs_Sz, mbs_r, mbs_sig);
	std::pair<std::vector<double>, std::vector<double>> Sta;
	std::pair <double, double> Ev;
	switch (method)
	{
	case 0: //Call Européen
		Ev = calcPayoffEUa();
		std::get<0>(Ev) = exp(-mbs_r*mbs_T)*std::get<0>(Ev);
		std::get<1>(Ev) = exp(-mbs_r*mbs_T)*std::get<1>(Ev);
		return Ev;
		break;
	case 1: //Call Asiatiaque
		std::get<0>(Sta).resize(Ti.size());
		std::get<1>(Sta).resize(Ti.size());
		Sta = ps1.CalculSta(Ti);
		std::get<0>(Ev) = exp(-mbs_r*mbs_T)*calcPayoffAS(std::get<0>(Sta));
		std::get<1>(Ev) = exp(-mbs_r*mbs_T)*calcPayoffAS(std::get<1>(Sta));
		return Ev;
		break;
	case 3: //Call down and in
		std::get<0>(Sta).resize(Ti.size());
		std::get<1>(Sta).resize(Ti.size());
		Sta = ps1.CalculSta(Ti);
		std::get<0>(Ev) = exp(-mbs_r*mbs_T)*calcPayoffUO(std::get<0>(Sta),M);
		std::get<1>(Ev) = exp(-mbs_r*mbs_T)*calcPayoffUO(std::get<1>(Sta),M);
		return Ev;
		break;
	default:
		return Ev;
	}
}

double MonteCarloBS::maxCondBrow(double ST)
{
	srand(time(NULL));
	double u = rand();
	return (ST + sqrt(ST*ST - 2*mbs_T*log(u)))/2;
}


double MonteCarloBS::calcC0(int P, std::vector<double> Ti) {
	PriceA ps1 = PriceA(mbs_Sz, mbs_r, mbs_sig);
	std::vector<double> St;
	St.resize(Ti.size());
	std::vector<double> H;
	H.resize(P);
	std::vector<double> Z;
	Z.resize(P);
	for (int p = 0; p < P; p++) {
		St = ps1.CalculSt(Ti);
		H[p] = exp(-mbs_r*mbs_T)*calcPayoffAS(St);
		Z[p] = St[Ti.size() - 1];
		//Z[p] = std::accumulate(St.begin(), St.end(), 0.0) / Ti.size();
	}
	double Hm = std::accumulate(H.begin(), H.end(), 0.0) / P;
	double Zm = std::accumulate(Z.begin(), Z.end(), 0.0) / P;
	double VZ = 0;

	double c = 0;
	for (int i = 0; i < P; i++) {
		c = c + (H[i] - Hm)*(Z[i] - Zm);
	}
	for (int j = 0; j < P; j++) {
		VZ = VZ + pow((Z[j] - Zm),2);
	}
	return c / VZ;
}


double MonteCarloBS::calcC1(int P, std::vector<double> Ti) {
	PriceA ps1 = PriceA(mbs_Sz, mbs_r, mbs_sig);
	std::vector<double> St;
	St.resize(Ti.size());
	std::vector<double> H;
	H.resize(P);
	std::vector<double> Z;
	Z.resize(P);
	for (int p = 0; p < P; p++) {
		St = ps1.CalculSt(Ti);
		H[p] = exp(-mbs_r*mbs_T)*calcPayoffAS(St);
		//Z[p] = St[Ti.size() - 1];
		Z[p] = std::accumulate(St.begin(), St.end(), 0.0) / Ti.size();
	}
	double Hm = std::accumulate(H.begin(), H.end(), 0.0) / P;
	double Zm = std::accumulate(Z.begin(), Z.end(), 0.0) / P;
	double VZ = 0;

	double c = 0;
	for (int i = 0; i < P; i++) {
		c = c + (H[i] - Hm)*(Z[i] - Zm);
	}
	for (int j = 0; j < P; j++) {
		VZ = VZ + pow((Z[j] - Zm), 2);
	}
	return c / VZ;
}

double MonteCarloBS::calcC2(int P, std::vector<double> Ti) {
	PriceA ps1 = PriceA(mbs_Sz, mbs_r, mbs_sig);
	std::vector<double> St;
	St.resize(Ti.size());
	std::vector<double> H;
	H.resize(P);
	std::vector<double> Z;
	Z.resize(P);
	for (int p = 0; p < P; p++) {
		St = ps1.CalculSt(Ti);
		H[p] = exp(-mbs_r*mbs_T)*calcPayoffAS(St);
		Z[p] = exp(-mbs_r*mbs_T)*std::max(0.0, St[Ti.size() - 1] - mbs_K);
	}
	double Hm = std::accumulate(H.begin(), H.end(), 0.0) / P;
	double Zm = std::accumulate(Z.begin(), Z.end(), 0.0) / P;
	double VZ = 0;

	double c = 0;
	for (int i = 0; i < P; i++) {
		c = c + (H[i] - Hm)*(Z[i] - Zm);
	}
	for (int j = 0; j < P; j++) {
		VZ = VZ + pow((Z[j] - Zm), 2);
	}
	return c / VZ;
}