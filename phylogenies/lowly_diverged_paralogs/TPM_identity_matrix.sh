#!/bin/bash
#SBATCH --mem=0
#SBATCH -t 72:00:00
#SBATCH -o sbatch_test.out
#SBATCH -p fn2

module purge
source /home/mjpete11/.bashrc
conda activate clustalo

# align SLC25 RNA sequences with muscle
muscle -in TPM.fa -clwstrictout TPM.clw

# generate identity matrix
clustalo -i TPM.clw --percent-id --distmat-out=TPM_pim.txt --full --force
