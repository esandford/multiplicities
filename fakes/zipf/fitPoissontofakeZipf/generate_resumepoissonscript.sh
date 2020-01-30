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
    
    scriptname=`echo "./resumefitPoissontofakeZipf_pop"$j"_run"$i".sh"`
    jobname=`echo "resumefitPoissontofakeZipf_pop"$j"_run"$i`

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
    echo "SBATCH --mem-per-cpu=3gb" >> $scriptname
    echo "" >> $scriptname
    echo "module load netcdf-gfortran/4.4.4" >> $scriptname
    echo "module load netcdf-fortran/4.4.4" >> $scriptname
    echo "" >> $scriptname
    echo "#./make.sh" >> $scriptname
    # for resume version
    nstepsalready=`wc -l "/rigel/home/es3197/multiPlanetDist_fixedForecaster/fakes/zipf/fitPoissontofakeZipf/mcmc_pop"$j"_run"$i".dat"`
    nstepsalready=`echo $nstepsalready | cut -d ' ' -f 1`
    echo "./resumesim "$i $nstepsalready $j >> $scriptname
    echo "" >> $scriptname
  done
done
echo 'Done generating sbatch scripts.'
