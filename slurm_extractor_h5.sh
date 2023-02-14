#!/bin/sh

# SLURM options:

#SBATCH --partition=htc               # Partition choice
#SBATCH --ntasks=1                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=9000                    # Memory in MB per default
#SBATCH --time=40:00:00

# Commands to be submitted:
module unload km3net_env
module load km3net_soft_env/1.9

python3 /sps/km3net/users/fbenfe/Moon_shadow/arca_h5_extractor.py ${FILE} ${OUTFILE}

#printenv | grep KM3
#python3 /sps/km3net/users/fbenfe/Moon_shadow/arca_data_extractor_fake_sources_single.py ${FILE} ${OUTFILE} > /dev/null
#echo "recreated cookie"
#python3 /sps/km3net/users/fbenfe/Moon_shadow/arca_data_extractor_fake_sources_single.py ${FILE} ${OUTFILE}

