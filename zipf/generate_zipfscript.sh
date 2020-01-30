#!/bin/bash

istart=0
iend=52

for i in `seq $istart $iend`; do
  echo ' '
  echo 'Generating sbatch script for run '$i' ('$(($i))' of '$(($iend))')'
  
  scriptname=`echo "./runZipf_"$i".sh"`
  jobname=`echo "runZipf_"$i`

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
  echo "./sim "$i >> $scriptname
  echo "" >> $scriptname

done
echo 'Done generating sbatch scripts.'
