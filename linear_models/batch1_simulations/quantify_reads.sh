#!/bin/bash
 
#SBATCH -p htc
#SBATCH -q normal
#SBATCH -o sbatch_test.out

module purge
source /home/mjpete11/.bashrc
conda activate data_science

Rscript num_reads_table.R
