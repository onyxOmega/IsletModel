#include "BetaCell.h"
#include "boost/random/random_device.hpp"
#include <random>
#include <ctime>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
//Randomizes cellular properties and saves to file
int main()
{
	srand (1);

//Initialize random variable matrix
double randomVars[cellNumber][12];
double cSum=0;

boost::mt19937 gen(4);
boost::random::uniform_real_distribution<> dis(0.25, 1.5);
//boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)

 //boost::normal_distribution<> gKATPv(2.57, 0.272);
 boost::normal_distribution<> gKATPv(2.31, .23);
boost::gamma_distribution<>gCoupS(4,4);


//Set coupling strengths here


//boost::normal_distribution<> gKtoarv(2.13, .232);

boost::normal_distribution<> gKtoarv(2.13, 0.231);
//boost::normal_distribution<> PCaERarv(0.096, .0096);
//boost::normal_distribution<> PCaERarv(0.15, 0.0156);

boost::normal_distribution<> PCaERarv(0.096, 0.0090);

//boost::normal_distribution<> pKATPv(1.0,0.000);


//boost::normal_distribution<> gGKCaBKv(2.3,0.2322);
boost::normal_distribution<> gGKCaBKv(2.31,0.230);

//boost::normal_distribution<> PNACAv(204,20);

boost::normal_distribution<> PNACAv(204,20);
//boost::normal_distribution<> Prelv(0.46,0.046);
 boost::normal_distribution<> Prelv(0.46,0.0446);

//boost::normal_distribution<> Popv(0.0005,0.00003);
 boost::normal_distribution<> Popv(0.0005,0.00000);
//boost::normal_distribution<> ATPv(4,0.4);
 boost::normal_distribution<> ATPv(4,0.4);
//boost::normal_distribution<> KBOXv(0.0000063,0.0000006);
 boost::normal_distribution<> KBOXv(0.0000063,0.0000006);

//boost::normal_distribution<> GLYCv(0.000126,0.0000315);
 boost::normal_distribution<> GLYCv(0.000126,0.0000315);
for(int i=0;i<cellNumber;i++)
{
//Calling the random number generators from BetaCell.h
randomVars[i][0]=gKATPv(gen);
randomVars[i][1]=gCoupS(gen);
cSum=cSum+randomVars[i][1];
randomVars[i][2]=gKtoarv(gen);
randomVars[i][3]=PCaERarv(gen);
randomVars[i][4]=gGKCaBKv(gen);
randomVars[i][5]=PNACAv(gen);
randomVars[i][6]=Prelv(gen);
randomVars[i][7]=Popv(gen);
randomVars[i][8]=ATPv(gen);
randomVars[i][9]=GLYCv(gen);
//Unique random number generator to help determine coupling conductance.
randomVars[i][10]=rand() % 1000000;
randomVars[i][11]=dis(gen);
}

double cMean=cSum/cellNumber;

////////////////Below sets the average coupling conductance///////////////
for(int i=0;i<cellNumber;i++)
{
randomVars[i][1]=randomVars[i][1]*0.00/cMean;
}
///////////////////////////////////////////////////////////////////
remove("RandomVars.txt");

//Write variables to file RandomVars.txt
ofstream outFile;
outFile.open("RandomVars.txt");
for(int i=0;i<cellNumber;i++)
{
for(int j=0;j<=11;j++)
{
outFile<<randomVars[i][j]<<' ';
}
outFile<<' '<<endl;
}




return 0;
}
