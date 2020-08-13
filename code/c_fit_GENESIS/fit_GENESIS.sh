#!/bin/bash -l

#SBATCH
#SBATCH --job-name=GENESIS
#SBATCH --time=72:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem 30G
# memory pool for all cores
#SBATCH --ntasks=1
# number of cpus (threads) per task (process)
#SBATCH --cpus-per-task=24
module load R/3.4.0
Rscript fit_GENESIS.R  $iter $pic
