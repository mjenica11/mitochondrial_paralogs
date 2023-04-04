#!/bin/bash

#SBATCH -o prob.out
#SBATCH -t 6-23:59:00 

module purge
source /home/mjpete11/.bashrc
conda activate peer
module load r/3.4.4

Rscript peer.R
