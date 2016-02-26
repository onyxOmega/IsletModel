Instructions:
Compiler essentials:

Linux: Needs build essential tools:

sudo apt-get install build-essential

Windows:

Cywgin (https://www.cygwin.com/) at minimum needs gcc-g++ package or GIT-Bash will work also. Can also be used on visual studio if you want to take that route (bigger pain)

Mac: Not 100% Sure

1) compile Sphere packing 

g++  neighbor.C spheres.C box.C sphere.C event.C heap.C read_input.C -o SpherePacking

2)Create a random Sphere Packing:

./SpherePacking sphereInput

3)Compile Islet generator

g++ GenerateSphere -o IsletGenerate

4)Create Islet of custom Size:

./IsletGenerate isletInput


5) compile Random Number Generator for cell heterogeneity

g++ -std=c++0x -I (direcotry to boost) -fopenmp RandomVars.cpp -o Random

6)Create Cellular Random Variables:

./Random 


7) Compile Simulation:

g++ -I (Directory to boost) MainFile.cpp -o Beta


8) Run Simulation

./Beta M0




Input files: 
sphereInput (tells sphere packing how many hard spheres to pack (set at 4000 currently)

isletInput (gives the name for matrix containing each cell's nearest neighboor and how many total cells will be used)

M0(input file for simulation containng output file names for potential, calcium etc.)

Boost (needed for ode solvers/random number generators etc..)
http://www.boost.org/

