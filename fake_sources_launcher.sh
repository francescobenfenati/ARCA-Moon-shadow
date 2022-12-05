#!/bin/bash
filename=$1
while read line; do
# reading each line
csv=$(echo $line | awk '{print $1}')
outfile_name=$(echo $line | awk '{print $2}')
beta_cut=$2
maxdeg=$3
bin_size=$4

job=$(echo $csv | rev | cut -d/ -f1 | rev)
job_name=$(echo ${job} | rev | cut -d. -f2 | rev)

#echo $csv
#echo $outfile_name
#echo $beta_cut
#echo $maxdeg
#echo $bin_size
#echo $job_name 

sbatch --job-name=${job_name} --output=/sps/km3net/users/fbenfe/Moon_shadow/logs/arca19/${job_name}.log --export=CSV=${csv},OUTFILE_NAME=/sps/km3net/users/fbenfe/Moon_shadow/combined_shadow/data/arca19/${outfile_name},BETA_CUT=${beta_cut},MAXDEG=${maxdeg},BIN_SIZE=${bin_size} slurm_analyzer.sh

done < $filename
