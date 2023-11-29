#!/bin/bash

#SBATCH -n 8
#SBATCH -o matrix.out

module purge
source /home/mjpete11/.bashrc
module load r/4.0.0

Rscript make_salmon_count_matrix.R
