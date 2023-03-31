#!/bin/bash
filename=$1
outfile=$2
logfile=$3

while read line; do
# reading each line
job=$(echo $line)
job_name=$(echo ${job} | rev | cut -d. -f2 | rev)

#echo $line
#echo $outfile
#echo $job_name 

sbatch --job-name=${job_name} --output=${logfile}/${job_name}.log --export=FILE=${line},OUTFILE=${outfile},KM3NET_DB_USERNAME=fbenfenati,KM3NET_DB_PASSWORD=Bull1sm0 /sps/km3net/users/fbenfe/Moon_shadow/slurm_extractor_fake_sources.sh


#OUTFILE_NAME=/sps/km3net/users/fbenfe/Moon_shadow/combined_shadow/${outfile_name},SIGMA=${sigma},BETA_CUT=${beta_cut},MAXDEG=${maxdeg},BIN_SIZE=${bin_size} slurm_analyzer.sh

done < $filename
