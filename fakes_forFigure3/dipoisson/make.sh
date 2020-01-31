#!/bin/bash
#

compiler='gfortran'
rm -f /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dipoisson/*.o /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dipoisson/*.mod /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dipoisson/sim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dipoisson/resumesim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dipoisson/*~ /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dipoisson/mcmc.dat
$compiler -O3 -c /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dipoisson/library.f90
$compiler -O3 -c /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dipoisson/model.f90
$compiler -O3 -o /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dipoisson/sim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/dipoisson/sim.f90 *.o

#$compiler -c -g -Wall -Wtabs -fcheck=all library.f90
#$compiler -c -g -Wall -Wtabs -fcheck=all model.f90
#$compiler -g -Wall -Wtabs -fcheck=all -o sim sim.f90 *.o
#$compiler -g -Wall -Wtabs -fcheck=all -o resumesim resumeSim.f90 *.o
