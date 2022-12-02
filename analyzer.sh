#!/bin/bash
filename=$1
while read line; do
# reading each line

job=$(echo $line | rev | cut -d/ -f1 | rev)
job_name=$(echo ${job} | rev | cut -d. -f2 | rev)
plot=$job_name"_2Dhist.pdf"
#echo $job
#echo $job_name
#echo $plot
sbatch --job-name=${job_name} --output=/sps/km3net/users/fbenfe/Moon_shadow/logs/${job_name}.log --export=FILE=${line},PLOT=/sps/km3net/users/fbenfe/Moon_shadow/combined_shadow/data/arca19/${plot} slurm_analyzer.sh

done < $filename
