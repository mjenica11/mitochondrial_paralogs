#!/bin/bash
#SBATCH -o PPTG.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate clustalo

# align SLC25 RNA sequences with muscle
#muscle -in PTTG.fa -clwstrictout PTTG.clw

# generate identity matrix
clustalo -i PTTG1.clw --percent-id --distmat-out=PTTG_pim.txt --full --force
