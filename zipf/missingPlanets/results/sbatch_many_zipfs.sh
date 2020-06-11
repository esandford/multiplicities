#!/bin/bash

./make.sh

istart=0
iend=10

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

  zipfscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/zipf/missingPlanets/runZipf_"$i".sh"`
  #zipfscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/zipf/missingPlanets/resumeZipf_"$i".sh"`
  
  sbatch $zipfscript

done
