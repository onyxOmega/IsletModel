#ifndef ISLETSIMULATOR_H
#define ISLETSIMULATOR_H

#include <string>

class IsletSimulator
{
	private:
	std::string userVarMatrix[2][10];
	double ktt, kdd, ktd, runTime, tStep;
	char const* timeOutput= "time.txt";						// time stamp for measuring progress.
	char const* potentialOutput= "potential.txt";			// membrane potential	
	char const* calciumOutput= "calcium.txt";			// intracellular calcium
	char const* sodiumOutput= "sodium.txt";				// intracellular sodium
	char const* potassiumOutput= "potassium.txt";	// intracellular potassium
	char const* caerOutput= "caer.txt";						// endoplasmic reticulum calcim
	char const* atpOutput= "atp.txt";						// intracellular ATP
	char const* adpOutput= "adp.txt";						// intracellular ADP
	char const* IRPOutput = "IRP.txt"; 					 	// immediately releasable pool		x[k+22] 
	char const* PPOutput = "PP.txt";							// primed pool								x[k+23] 
	char const* DPOutput ="DP.txt";							// docked pool								x[k+24]
	char const* FIPOutput = "FIP.txt";						// fused pool (FHP in the paper)
	char const* RIPOutput = "RIP.txt";						// releasing pool (RHP in the paper)
	char const* capOutput= "cap.txt";						// ??? I don't know what this one is, it's Variable 28 in the X vector
	char const* noiseOutput= "noise.txt";					// ??? I don't know how this works, var 29 in X vector
	char const* varFileName = "UserDefinedVars.txt";
		
	public:
		void initialize(const char *);
		void setUserDefinedVars();
		void setDefaultVars();
		void purgeOutputs();
		double get_ktt() const;
		double get_kdd() const;
		double get_ktd() const;
		double get_runTime() const;
};
#endif