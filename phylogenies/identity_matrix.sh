#!/bin/bash
#SBATCH --mem=0
#SBATCH -t 72:00:00
#SBATCH -o sbatch_test.out
#SBATCH -p fn2

module purge
source /home/mjpete11/.bashrc
conda activate clustalo

# align SLC25 DNA sequences with muscle
muscle -in HGNC_named_combined.fasta -clwstrictout aligned.clw

# generate identity matrix
clustalo -i aligned.clw --percent-id --distmat-out=pim.txt --full --force
