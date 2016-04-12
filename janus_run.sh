#!/bin/bash

#Setting the name of the job
#PBS -N IsletSimObjectiveTest
#Setting a walltime for the job
#PBS -l walltime=1:00:00
#Selecting processors
#PBS -l nodes=1:ppn=12
#SBATCH --reservation=janus-serial
 
cd $PBS_O_WORKDIR

#Execute

./Beta.exe default >> runtimeOutput.txt