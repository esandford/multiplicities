#!/bin/bash

./make.sh

istart=0
iend=52

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

  dipoissonscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/dipoisson/runDipoisson_"$i".sh"`
  #dipoissonscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/dipoisson/resumeDipoisson_"$i".sh"`
  
  sbatch $dipoissonscript

done
