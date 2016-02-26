#include "Islet.h"
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
	IsletSimulator::setUserDefinedVars();
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
	}
	
	/* Set default values for anything that wasn't initialized in the
		previous loop. -WLF
	*/
	if(!ktt)
	{
		ktt = 0.05;
	}
	if(!kdd)
	{
		kdd = 0.01;
	}
	if(!ktd)
	{
		ktd = 0.026;
	}
	cout << "ktt = " << ktt << "; kdd = " << kdd << "; ktd = " << ktd;
}

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
