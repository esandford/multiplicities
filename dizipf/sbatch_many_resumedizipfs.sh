#!/bin/bash

./make.sh

istart=0
iend=52

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

  #dizipfscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/dizipf/runDizipf_"$i".sh"`
  dizipfscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/dizipf/resumeDizipf_"$i".sh"`
  
  sbatch $dizipfscript

done
