#!/bin/bash


istart=0
iend=99

model="zipf"

makescript=`echo "./"$model"/make.sh"`
simscript=`echo "./"$model"/sim"`
chain=`echo "../"$model"Chain.dat"`

. "$makescript"

for i in `seq $istart $iend`; do

    #get a random index between 1 and 50001, inclusive 
    randIdx=`awk 'BEGIN {
       # seed
       srand()
       print int(1+rand() * 50001)
    }'`

    #echo $randIdx
    echo $i

    radfile=`echo "/Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/"$model"/fake_keplerradii_"$i".dat"`
    detfile=`echo "/Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/"$model"/fake_keplerdetections_"$i".dat"`
    countfile=`echo "/Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/"$model"/fake_keplercounts_"$i".dat"`
    
    $simscript $chain $randIdx $radfile $detfile $countfile

done