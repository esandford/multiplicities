#!/bin/bash

./make.sh

istart=0
iend=52

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

  diuniformscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/diuniform/runDiuniform_"$i".sh"`
  #diuniformscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/diuniform/resumeDiuniform_"$i".sh"`
  
  sbatch $diuniformscript

done
