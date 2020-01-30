#!/bin/bash

./make.sh

istart=0
iend=52

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

  #constantscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/constant/runConstant_"$i".sh"`
  constantscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/constant/resumeConstant_"$i".sh"`
  
  sbatch $constantscript

done
