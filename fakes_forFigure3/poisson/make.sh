#!/bin/bash
#

compiler='gfortran'
rm -f /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/poisson/*.o /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/poisson/*.mod /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/poisson/sim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/poisson/resumesim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/poisson/*~ /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/poisson/mcmc.dat
$compiler -O3 -c /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/poisson/library.f90
$compiler -O3 -c /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/poisson/model.f90
$compiler -O3 -o /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/poisson/sim /Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/fakes_forFigure3/poisson/sim.f90 *.o

#$compiler -c -g -Wall -Wtabs -fcheck=all library.f90
#$compiler -c -g -Wall -Wtabs -fcheck=all model.f90
#$compiler -g -Wall -Wtabs -fcheck=all -o sim sim.f90 *.o
#$compiler -g -Wall -Wtabs -fcheck=all -o resumesim resumeSim.f90 *.o
