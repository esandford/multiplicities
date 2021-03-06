#!/bin/bash

# number of mcmc runs
istart=0
iend=10

# number of fake kepler populations
jstart=0
jend=9

for i in `seq $istart $iend`; do
  for j in `seq $jstart $jend`; do
    echo ' '
    echo 'Generating sbatch script for run '$i' ('$(($i))' of '$(($iend))')'
    
    scriptname=`echo "./runfitPoissontofakePoisson_pop"$j"_run"$i".sh"`
    jobname=`echo "runfitPoissontofakePoisson_pop"$j"_run"$i`

    # Generate empty file to hold sbatch script
    rm $scriptname
    touch $scriptname

    echo "#!/bin/sh" >> $scriptname
    echo "#" >> $scriptname
    echo "#" >> $scriptname
    echo "#SBATCH --account=astro" >> $scriptname
    echo "#SBATCH --job-name="$jobname >> $scriptname
    echo "#SBATCH -c 1" >> $scriptname
    echo "#SBATCH --time=120:00:00" >> $scriptname
    echo "#SBATCH --mem-per-cpu=3gb" >> $scriptname
    echo "" >> $scriptname
    echo "module load netcdf-gfortran/4.4.4" >> $scriptname
    echo "module load netcdf-fortran/4.4.4" >> $scriptname
    echo "" >> $scriptname
    echo "#./make.sh" >> $scriptname
    echo "./sim "$i" "$j >> $scriptname
    echo "" >> $scriptname
  done
done
echo 'Done generating sbatch scripts.'
