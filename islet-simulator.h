/* Constructor for IsletSimulator custom class.

	This class stores islet specific data (including an array of Beta
	Cell objects within the simulated islet), and calculates property
	changes over time using a euler first order linar approximation on
	a series of interdependant ODEs.
	
	Authors: William Fischer, Matt Wescott
*/

#ifndef ISLETSIMULATOR_H
#define ISLETSIMULATOR_H

#include "beta-cell.h"

#include <string>
#include <vector>

using namespace std;

class IsletSimulator
{
	private:
		// User defined variables
		std::string userVarMatrix[2][10];
		double ktt, kdd, ktd, runTime, stepTime;
		
		/* Built in universal variables (for whole islet or uniform across
		    all beta cells)
		*/
		
		// Beta Cell vector
		vector<BetaCell> betaCellVector;
		
	public:
		IsletSimulator();										// Class constructor
		void initialize(const char *);
		void setDefaultVars();
		void setUserDefinedVars();
		double get_ktt() const;
		double get_kdd() const;
		double get_ktd() const;
		double get_runTime() const;
		double get_stepTime() const;
};

#endif