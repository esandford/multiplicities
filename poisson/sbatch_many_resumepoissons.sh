#!/bin/bash

./make.sh

istart=0
iend=52

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

  #poissonscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/poisson/runPoisson_"$i".sh"`
  poissonscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/poisson/resumePoisson_"$i".sh"`
  
  sbatch $poissonscript

done
