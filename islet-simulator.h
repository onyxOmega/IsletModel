/* Constructor for IsletSimulator custom class.

	This class stores islet specific data (including an array of Beta
	Cell objects within the simulated islet), and calculates property
	changes over time using a euler first order linar approximation on
	a series of interdependant ODEs.
	
	Authors: William Fischer, Matt Wescott
*/

#ifndef ISLETSIMULATOR_H
#define ISLETSIMULATOR_H

#include "islet-file-handler.h"
#include "islet-data-structures.h"


#include <string>
#include <vector>


using namespace std;

class IsletSimulatorClass
{
	private:
		IsletFileHandlerClass fileHandler;
		
		// User defined variables
		std::string userVarMatrix[2][10];
		double runTime, stepTime;
		
		// other variables
		double Glucose;
		int cellNumber;
		double Icoup;
		
		// Beta Cell vector
		IsletStructure islet;
		vector<BetaCellStructure> betaCells;
		
		const double R = 8.3143;
		const double Tem = 310.15;
		const double F = 96.4867;
		const double RTF = R*Tem/F;
		const double RTF2 = RTF/2;
		const double yini0 = -69.8663703359279;
		const double yini1 = 7.92502913389466;
		const double yini2 = 125.248586232226;
		const double yini3 = 7.56138347594955E-05;
		const double yini4 = 0.0047417301517704;

		
	public:
		IsletSimulatorClass(IsletFileHandlerClass);					// Class constructor
		void setDefaultVars();
		void setInitialBetaCellVars();
		void setUserDefinedVars();
		void simulationLoop();
		void setNearestNeighbors();
		double get_ktt() const;
		double get_kdd() const;
		double get_ktd() const;
		double get_runTime() const;
		double get_stepTime() const;
};

#endif