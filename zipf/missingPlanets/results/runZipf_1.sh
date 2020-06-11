#!/bin/sh
#
#
#SBATCH --account=astro
#SBATCH --job-name=runZipf_1
#SBATCH -c 1
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=3gb

module load netcdf-gfortran/4.4.4
module load netcdf-fortran/4.4.4

#./make.sh
./sim 1

