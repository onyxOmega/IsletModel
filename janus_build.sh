#!/bin/bash
module load intel
icpc -std=c++0x -I /projects/fischwil/boost_1_54_0/ -fopenmp main.cpp islet-simulator.cpp islet-file-handler.cpp ch-r2.cpp -o simulate.exe