/* Custom Class: IsletFileHandler

	Handles input and output files used in the Islet of Langerhans 
	computer model.
	
	Authors: William Fischer, Matt Wescott
*/

#ifndef ISLETFILEHANDLER_H
#define ISLETFILEHANDLER_H

#include <vector>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/detail/config.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include <boost/math/distributions/skew_normal.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;

class IsletFileHandler 
{
	private:
		// Output filenames.
		char const* timeOutput = "output/time.txt";				
		char const* potentialOutput = "output/potential.txt";			// membrane potential	
		char const* calciumOutput = "output/calcium.txt";			// intracellular calcium
		char const* sodiumOutput = "output/sodium.txt";				// intracellular sodium
		char const* potassiumOutput = "output/potassium.txt";	// intracellular potassium
		char const* caerOutput = "output/caer.txt";						// endoplasmic reticulum calcim
		char const* atpOutput = "output/atp.txt";						// intracellular ATP
		char const* adpOutput = "output/adp.txt";						// intracellular ADP
		char const* IRPOutput = "output/IRP.txt"; 					 	// immediately releasable pool		x[k+22] 
		char const* PPOutput = "output/PP.txt";							// primed pool								x[k+23] 
		char const* DPOutput = "output/DP.txt";							// docked pool								x[k+24]
		char const* FIPOutput = "output/FIP.txt";						// fused pool (FHP in the paper)
		char const* RIPOutput = "output/RIP.txt";						// releasing pool (RHP in the paper)
		char const* capOutput = "output/cap.txt";						// ??? I don't know what this one is, it's Variable 28 in the X vector
		char const* noiseOutput = "output/noise.txt";					// ??? I don't know how this works, var 29 in X vector
		
		// Input filenames. Runtime modifications not yet implemented
		char const* userVarsFile = "input/UserDefinedVars.txt";
		char const* cellPropertiesFile = "input/vars5exo.txt";
		char const* nnFile = "input/NN10A.txt";
		char const* randomVarsFile = "input/RandomVars.txt";
		char const* cellPostitionFile = "input/XYZpos.txt";
	
	public:
		void initialize(const char*);
		void writeOutputs(vector_type, int);
		
		void purgeOutputFiles();
		char const* get_userVarsFile();
		char const* get_cellPropertiesFile();
		char const* get_cellPositionFile();
		char const* get_nnFile();
		char const* get_randomVarsFile();
};

# endif