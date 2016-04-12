#!/bin/bash

g++ -std=c++0x -I /home/boost_1_54_0/ -fopenmp main.cpp islet-simulator.cpp islet-file-handler.cpp ch-r2.cpp -o Beta.exe
./Beta.exe input/UserDefinedVars.txt
