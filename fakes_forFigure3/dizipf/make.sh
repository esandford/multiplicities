#!/bin/bash
#

compiler='gfortran'
rm -f /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dizipf/*.o /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dizipf/*.mod /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dizipf/sim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dizipf/resumesim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dizipf/*~ /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dizipf/mcmc.dat
$compiler -O3 -c /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dizipf/library.f90
$compiler -O3 -c /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dizipf/model.f90
$compiler -O3 -o /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dizipf/sim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dizipf/sim.f90 *.o

#$compiler -c -g -Wall -Wtabs -fcheck=all library.f90
#$compiler -c -g -Wall -Wtabs -fcheck=all model.f90
#$compiler -g -Wall -Wtabs -fcheck=all -o sim sim.f90 *.o
#$compiler -g -Wall -Wtabs -fcheck=all -o resumesim resumeSim.f90 *.o
