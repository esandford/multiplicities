#!/bin/bash

./make.sh

istart=0
iend=52

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Submitting run '$i' ('$(($i))' of '$(($iend))')'

  #diconstantscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/diconstant/runDiconsant_"$i".sh"`
  diconstantscript=`echo "/rigel/home/es3197/multiPlanetDist_fixedForecaster/diconstant/resumeDiconstant_"$i".sh"`
  
  sbatch $diconstantscript

done
