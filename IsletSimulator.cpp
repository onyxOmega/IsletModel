#include "IsletSimulator.h"
#include <iostream>
#include <string>
#include <stdio.h>
#include <ostream>
#include <fstream>
#include <sstream>

#include <boost/lexical_cast.hpp>

using namespace std;

/* This is an incomplete initialization of the Islet. I plan to gradually
	implement functionality, transferring one thing at a time from the 
	original code base to objective code and testing outputs. -WLF
*/

void IsletSimulator::initialize(const char * fileChars)
{	
	/* First part of initialization, take string literal output filenames, add
		a user defined path, and generate char pointers from them. 
		-WLF 
	*/

	/* The following is a general implementation for pulling in a series 
		of user defined variable values from an imput file. Not yet 
		implemented for variable assignment. So it doesn't do anything
		useful yet. - WLF
 	*/
	
	ifstream userVarFile;
	userVarFile.open(fileChars);
	
	for (int i = 0; !userVarFile.eof(); i++)
	{	
		/* The user defined variables are passed into a variable matrix
			"userVarMatrix" from the input file as string types. They are 
			error checked, converted into the appropriate data type, and 
			assigned in "setUserDefinedVars()". -WLF
		*/
		char buffer[20];
		stringstream varStream;
		userVarFile.getline(buffer, 20);
		varStream << buffer;
		varStream.getline(buffer, 10, '=');
		string strBuffer(buffer);
		userVarMatrix[0][i] = strBuffer;
		varStream.getline(buffer, 10, ';');
		strBuffer.assign(buffer);
		userVarMatrix[1][i] = strBuffer;
	}
	IsletSimulator::setDefaultVars();
	IsletSimulator::setUserDefinedVars();
}

void IsletSimulator::setDefaultVars()
{
	runTime = 500.0;
	ktt = 0.05;
	kdd = 0.01;
	ktd = 0.026;
}
/*	This function sets user defined variable values for any variables
listed in the imput file, and sets default values for any variable not
listed. Allows control of runtime parameters without altering the code
each run.
Currently implemented for ktt, ktd, kdd.
-WLF
*/

void IsletSimulator::setUserDefinedVars()
{

	for(int i = 0; userVarMatrix[0][i] != ""; i++)
	{
		/* If a given variable is on the list, set the value. Used boost
			lexical_cast to convert string to double because stod, atof,
			and strod aren't correctly implemented in cygwin. -WLF
		*/
		if(userVarMatrix[0][i] == "ktt")
		{
			ktt = boost::lexical_cast<double>(userVarMatrix[1][i]);			
		}
		if(userVarMatrix[0][i] == "kdd")
		{
			kdd = boost::lexical_cast<double>(userVarMatrix[1][i]);			
		}
		if(userVarMatrix[0][i] == "ktd")
		{
			ktd = boost::lexical_cast<double>(userVarMatrix[1][i]);			
		}
		if(userVarMatrix[0][i] == "runTime")
		{
			runTime = boost::lexical_cast<double>(userVarMatrix[1][i]);			
		}
	}
	
	/* Set default values for anything that wasn't initialized in the
		previous loop. -WLF
	*/
}

// Getters
double IsletSimulator::get_ktt() const
{
	return ktt;
}

double IsletSimulator::get_kdd() const
{
	return kdd;
}

double IsletSimulator::get_ktd() const
{
	return ktd;
}

double IsletSimulator::get_runTime() const
{
	return runTime;
}


//	delete all the output files in a given path
void IsletSimulator::purgeOutputs()
{
	if (remove(potentialOutput)) perror("Error 5");
	if (remove(timeOutput)) perror("Error 6");
	if (remove(calciumOutput)) perror("Error 7");
	if (remove(sodiumOutput)) perror("Error 8");
	if (remove(potassiumOutput)) perror("Error 9");
	if (remove(caerOutput)) perror("Error 10");
	if (remove(atpOutput)) perror("Error 11");
	if (remove(adpOutput)) perror("Error 12");
	if (remove(PPOutput)) perror("Error 13");
	if (remove(IRPOutput)) perror("Error 14");
	if (remove(DPOutput)) perror("Error 15");
	if (remove(FIPOutput)) perror("Error 16");
	if (remove(RIPOutput)) perror("Error 17");
	if (remove(capOutput)) perror("Error 18");
	if (remove(noiseOutput)) perror("Error 19");
}
