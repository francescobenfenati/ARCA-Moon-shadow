#!/bin/sh

# SLURM options:

#SBATCH --partition=htc               # Partition choice
#SBATCH --ntasks=1                    # Run a single task (by default tasks == CPU)
#SBATCH --mem=9000                    # Memory in MB per default


# Commands to be submitted:
module unload km3net_env
module load km3net_soft_env/1.9

./contour_data.exe ${CSV} ${OUTFILE_NAME} ${BETA_CUT} ${MAXDEG} ${BIN_SIZE}
