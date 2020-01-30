#!/bin/bash

./make.sh

istart=0
iend=52

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

  uniformscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/uniform/runUniform_"$i".sh"`
  #uniformscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/uniform/resumeUniform_"$i".sh"`
  
  sbatch $uniformscript

done
