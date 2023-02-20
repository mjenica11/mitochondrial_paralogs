#!/bin/bash

#SBATCH -n 8
#SBATCH -o sbatch_combat_seq.out
#SBATCH -p fn3 

module purge
module load r/4.0.0

Rscript combat_seq.R
