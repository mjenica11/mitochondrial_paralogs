#!/bin/bash
#SBATCH --mem=0
#SBATCH -o pims.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate clustalo

# generate identity matrix
clustalo -i GMCL1.clw --percent-id --distmat-out=GMCL1_pim.txt --full --force
clustalo -i PTTG1.clw --percent-id --distmat-out=PTTG1_pim.txt --full --force
clustalo -i TPM.clw --percent-id --distmat-out=TPM_pim.txt --full --force
