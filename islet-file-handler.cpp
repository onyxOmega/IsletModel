#include "islet-file-handler.h"
#include <iostream>
#include <string>
#include <stdio.h>
#include <ostream>
#include <fstream>
#include <sstream>
#include <fstream>
#include <vector>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;
using namespace std;

void IsletFileHandler::initialize(char const* outputPath)
{
	
}

char const* IsletFileHandler::get_userVarsFile()
{
	return userVarsFile;
}

char const* IsletFileHandler::get_cellPropertiesFile()
{
	return cellPropertiesFile;
}

char const* IsletFileHandler::get_nnFile()
{
	return nnFile;
}

char const* IsletFileHandler::get_randomVarsFile()
{
	return randomVarsFile;
}

void IsletFileHandler::writeOutputs(vector_type x, int cellNumber)
{
	ofstream outfilePotential;
	outfilePotential.open(potentialOutput,ios::app);
	ofstream outfileCalcium;
	outfileCalcium.open(calciumOutput,ios::app);
	ofstream outfileSodium;
	outfileSodium.open(sodiumOutput,ios::app);
	ofstream outfilePotassium;
	outfilePotassium.open(potassiumOutput,ios::app);
	ofstream outfileCaer;
	outfileCaer.open(caerOutput,ios::app);
	ofstream outfileATP;
	outfileATP.open(atpOutput,ios::app);
	ofstream outfileADP;
	outfileADP.open(adpOutput,ios::app);
	//ofstream outfileO1;
	//ofstream outfileO2;
	//ofstream outfileC1;
	//ofstream outfileC2;
	//outfileO1.open(O1Output,ios::app);
	//outfileO2.open(O2Output,ios::app);
	//outfileC1.open(C1Output,ios::app);
	//outfileC2.open(C2Output,ios::app);
	ofstream outfileIRP;
	outfileIRP.open(IRPOutput,ios::app);
	ofstream outfilePP;
	outfilePP.open(PPOutput,ios::app);
	ofstream outfileDP;
	outfileDP.open(DPOutput,ios::app);
	ofstream outfileFIP;
	outfileFIP.open(FIPOutput,ios::app);
	ofstream outfileRIP;
	outfileRIP.open(RIPOutput,ios::app);
	ofstream outfileCap;
	outfileCap.open(capOutput,ios::app);
	ofstream outfileNoise;
	outfileNoise.open(noiseOutput,ios::app);
	
	for (int k=0;k<30*cellNumber;k=k+30)
	{
		outfilePotential<<x[k]<<' ';
		outfileSodium<<x[k+1]<<' ';
		outfilePotassium<<x[k+2]<<' ';
		outfileCaer<<x[k+4]<<' ';
		outfileCalcium<<x[k+3]<<' ';
		outfileATP<<x[k+5]<<' ';
		outfileADP<<x[k+6]<<' ';
		//outfileO1<<x[k+18]<<' ';
		//outfileO2<<x[k+19]<<' ';
		//outfileC1<<x[k+20]<<' ';
		//outfileC2<<x[k+21]<<' ';
		outfileIRP<<x[k+22]<<' ';
		outfilePP<<x[k+23]<<' ';
		outfileDP<<x[k+24]<<' ';
		outfileFIP<<x[k+26]<<' ';
		outfileRIP<<x[k+27]<<' ';
		outfileCap<<x[k+28]<<' ';
		outfileNoise<<x[k+29]<<' ';
	}
	
	outfilePotential<<' '<<endl;
	outfileCalcium<<' '<<endl;
	outfileSodium<<' '<<endl;
	outfilePotassium<<' '<<endl;
	outfileCaer<<' '<<endl;
	outfileATP<<' '<<endl;
	outfileADP<<' '<<endl;
	//outfileO1<<' '<<endl;
	//outfileO2<<' '<<endl;
	//outfileC1<<' '<<endl;
	//outfileC2<<' '<<endl;
	outfileIRP<<' '<<endl;
	outfilePP<<' '<<endl;
	outfileDP<<' '<<endl;
	outfileFIP<<' '<<endl;
	outfileRIP<<' '<<endl;
	outfileCap<<' '<<endl;
	outfileNoise<<' '<<endl;
	outfilePotential.close();
	outfileATP.close();
	outfileADP.close();
	//outfileO1.close();
	//outfileO2.close();
	//outfileC1.close();
	//outfileC2.close();
	outfileCalcium.close();
	outfileSodium.close();
	outfilePotassium.close();
	outfileCaer.close();
	outfileIRP.close();
	outfilePP.close();
	outfileDP.close();
	outfileFIP.close();
	outfileRIP.close();
	outfileCap.close();
	outfileNoise.close();
}
			
void IsletFileHandler::purgeOutputFiles()
{
	if (remove(potentialOutput)) perror("Error 5");
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