/*	Custom data structures "IsletStructure" and "BetaCellStructure"
		for use in Islet simulation model.
	
		Author: William Fischer
*/

#ifndef ISLETDATASTRUCTURES_H
#define ISLETDATASTRUCTURES_H

#include <string>
#include <vector>

using namespace std;

struct IsletStructure
{	
	/* Data structure that stores islet parameters, or Beta Cell
		parameters that are constant across all beta cells. Items can
		be moved from this structure to the Beta Cell structure in order
		to use different values for different beta cells. This structure was
		made to avoid re-defining things for every cell at every time
		step. -WLF
	*/
	
	double ktt, kdd, ktd;
	
	double KRe;
	double Kfa;
	double Stoichi;
	double Rvol;
	double kATPCa;
	double kATP;
	double kADPf;
	double kADPb;
	double Naout;
	double Kout;
	double Caout;
	double Cm;
	double volER;
	double fi;
	double fer;
	double RCaLNa;
	double RCaLK;
	double PCaL;
	
	// Ikslow
	double PKslow;
	double nKslow;
	double KdKslow;

	// Ikdr
	double pKDr;
	double P_PMCA;
	double K_PMCA;

	// ICRAN params
	double PCRAN;
	double KCaer;
	double RNa_K_CRAN;
	double pIbNSC;
	double KTRPM;
	double pTRPM;
	double RNa_K_TRPM;

	//Ikatp params
	double gKATP;
	
	// NA/Ca exchange params
	double KdNao;
	double KdCao;
	double KdNai;
	double KdCai;
	double k3;
	double k4;
	
	// Na/K pump
	double Pii;
	double Proton;
	double Kd_MgATP;
	double Kd_Nao0;
	double Kd_Nai0;
	double Kd_Ko0;
	double Kd_Ki0;
	double delta_Nao;
	double delta_Nai;
	double delta_Ko;
	double delta_Ki;
	double k1_plus;
	double k2_plus;
	double k3_plus;
	double k4_plus;
	double k1_minus;
	double k2_minus;
	double k3_minus;
	double k4_minus;
	double PNaK;
	
	//Metabolism
	double Nt;
	double KCarp;
	
	//Glycolysis and oxidative phosph;
	double KmATP;
	double hgl;
	double Kg;
	
	//Check the Pop value
	double Kop;
	
	// IKATP: ATP gated potassium channel current:
	double residual;
	double kPrime;
	
	double taup;
	
	long double fusionMax;
	long double nFuse;
	long double K_I;
	
	// rate into/out of immediately releasable primed and docked 
	// pools, and from reserve
	long double r1;
	long double r_1;
	long double r2;
	long double r_2;
	long double Rres;
	long double R_res;
	
	// for cAMP later
	long double CaN;
	long double Kp;
	long double cN;
	long double cKi;
	
	//rate of movement from fusion to release pool (ms)
	long double u2;
	
	//rate of release (ms)
	long double u3;
	long double F_md;	
};

struct BetaCellStructure
{	
	/* Data structure that stores islet parameters, or Beta Cell
		parameters that are constant across all beta cells. -WLF
	*/
	
	/*These arrays are used to simplify the code block that calculates
		dxdt and the linear approximation steps
	*/
	double x[30];
	double dxdt[30];
		
	// Array of nearest neighbor cell numbers	
	vector<int> nnVector;				// list of nearest neighbor (coupled) cell index values
	int nnCount;								// total number of coupled cells.
	
	// variables populated from the randomVarsFile, in order of population
	double gKATPar;
	double gCoup;
	double gKtoar;
	double PCaER;
	double gKCaBKar;
	double PNACAar;
	double Prelar;
	double Popar;
	double ATPar;
	double KRev;
	double RandomSeed;
	double gChR2;
	
	// Variables populated from the varsFile, in order of population.
	/*
	double Vm;																				// x[0]
	double Nai;																				// x[1]
	double Ki;																				// x[2]
	double Cai;																				// x[3]
	double Caer;																			// x[4]
	double ATP;																			// x[5]
	double MgADP;																		// x[6]
	double Re;																				// x[7]
	double q_KDr;																		// x[8]
	double d_CaL;																		// x[9]
	double U_CaL;																		// x[10]
	double fus;																				// x[11]
	double p_KDr;																		// x[12]
	double m_Kto;																		// x[13]
	double h_Kto;																			// x[14]
	double E1_tota;																		// x[15]
	double I1;																				// x[16]
	double I2;																				// x[17]
	double O1;																				// x[18]
	double O2;																				// x[19]
	double C1;																				// x[20]
	double C2;																				// x[21]
	double IRP;																			// x[22]
	double PP;																				// x[23]
	double DP;																				// x[24]
	double RES;																			// x[25]
	double FIP;																			// x[26]
	double RIP;																			// x[27]
	double Cap;																			// x[28]
	double Pns;																			// x[29]
	*/
	
	// Other
	double IKATPvec;
};

#endif