#!/bin/bash

./make.sh

istart=0
iend=52

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

  diexponentialscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/diexponential/runDiexponential_"$i".sh"`
  #diexponentialscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/diexponential/resumeDiexponential_"$i".sh"`
  
  sbatch $diexponentialscript

done
