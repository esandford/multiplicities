#!/bin/bash
#

compiler='gfortran'
rm -f /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/constant/*.o /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/constant/*.mod /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/constant/sim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/constant/resumesim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDis_fixedForecaster/fakes_forFigure3/constant/*~ /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/constant/mcmc.dat
$compiler -O3 -c /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/constant/library.f90
$compiler -O3 -c /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/constant/model.f90
$compiler -O3 -o /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/constant/sim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/constant/sim.f90 *.o

#$compiler -c -g -Wall -Wtabs -fcheck=all library.f90
#$compiler -c -g -Wall -Wtabs -fcheck=all model.f90
#$compiler -g -Wall -Wtabs -fcheck=all -o sim sim.f90 *.o
#$compiler -g -Wall -Wtabs -fcheck=all -o resumesim resumeSim.f90 *.o
