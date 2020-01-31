#!/bin/bash
#

compiler='gfortran'
rm -f /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/diuniform/*.o /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/diuniform/*.mod /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/diuniform/sim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/diuniform/resumesim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/diuniform/*~ /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/diuniform/mcmc.dat
$compiler -O3 -c /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/diuniform/library.f90
$compiler -O3 -c /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/diuniform/model.f90
$compiler -O3 -o /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/diuniform/sim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/diuniform/sim.f90 *.o

#$compiler -c -g -Wall -Wtabs -fcheck=all library.f90
#$compiler -c -g -Wall -Wtabs -fcheck=all model.f90
#$compiler -g -Wall -Wtabs -fcheck=all -o sim sim.f90 *.o
#$compiler -g -Wall -Wtabs -fcheck=all -o resumesim resumeSim.f90 *.o
