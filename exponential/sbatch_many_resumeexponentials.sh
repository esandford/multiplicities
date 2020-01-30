#!/bin/bash

./make.sh

istart=0
iend=52

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

  #exponentialscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/exponential/runExponential_"$i".sh"`
  exponentialscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/exponential/resumeExponential_"$i".sh"`
  
  sbatch $exponentialscript

done
