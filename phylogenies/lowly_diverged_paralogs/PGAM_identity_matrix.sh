#!/bin/bash
#SBATCH --mem=0
#SBATCH -o pgam.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate clustalo

# align SLC25 RNA sequences with muscle
muscle -in PGAM.fa -clwstrictout PGAM.clw

# generate identity matrix
clustalo -i PGAM.clw --percent-id --distmat-out=PGAM_pim.txt --full --force
