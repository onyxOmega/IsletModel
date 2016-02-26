#!/bin/bash

#Setting the name of the job

#PBS -N Beta_Cell_Cluster_Ramp

#Setting a walltime for the job

#PBS -l walltime=7:00:00

#Selecting processors

#PBS -l nodes=1:ppn=12



#SBATCH --reservation=janus-serial 
cd $PBS_O_WORKDIR


#Execute

g++ -std=c++0x -I /projects/notarya/boost_1_54_0/ -fopenmp RandomVars.cpp -o RandomGenerator
./RandomGenerator

g++ -I /projects/notarya/boost_1_54_0/ -fopenmp MainFile.cpp -o Beta
./Beta 0B >> time_output0.txt

