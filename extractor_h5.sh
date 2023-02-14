#!/bin/bash
filename=$1
outfile=$2
logfile=$3
while read line; do
# reading each line

job=$(echo $line)
job_name=$(echo ${job} | rev | cut -d_ -f1 | rev | cut -d. -f1)

#echo $job
#echo $job_name

sbatch --job-name=${job_name} --output=${logfile}/${job_name}.log --export=FILE=${line},OUTFILE=${outfile},KM3NET_DB_USERNAME=fbenfenati,KM3NET_DB_PASSWORD=Bull1sm0 /sps/km3net/users/fbenfe/Moon_shadow/slurm_extractor_h5.sh
done < $filename
