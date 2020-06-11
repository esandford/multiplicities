#!/bin/bash
#

compiler='gfortran'
rm -f *.o *.mod sim resumesim *~ mcmc.dat
#$compiler -O3 -c library.f90
#$compiler -O3 -c model.f90
#$compiler -O3 -o sim sim.f90 *.o

$compiler -c -g -Wall -Wtabs -fcheck=all library.f90
$compiler -c -g -Wall -Wtabs -fcheck=all model.f90
$compiler -g -Wall -Wtabs -fcheck=all -o sim sim.f90 *.o
#$compiler -g -Wall -Wtabs -fcheck=all -o resumesim resumeSim.f90 *.o
