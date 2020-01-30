#!/bin/bash

./make.sh

# number of mcmc runs
istart=0
iend=10

# number of fake kepler populations
jstart=0
jend=9

for i in `seq $istart $iend`; do
  for j in `seq $jstart $jend`; do
	echo ' '
	echo ' '
	echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

	poissonscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/fakes/poisson/fitPoissontofakePoisson/resumefitPoissontofakePoisson_pop"$j"_run"$i".sh"`
	sbatch $poissonscript

  done
done
