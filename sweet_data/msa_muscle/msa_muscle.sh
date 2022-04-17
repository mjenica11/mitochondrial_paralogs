#!/bin/bash
#SBATCH -o sbatch.out

module purge
source /home/mjpete11/.bashrc
conda activate clustalo

# align unmapped reads of a non-failing heart control with muscle
muscle -in unmapped_non_failing.fa -clwstrictout unmapped_non_failing.clw

# align mapped reads of a non-failing heart control with muscle
muscle -in mapped_non_failing.fa -clwstrictout mapped_non_failing.clw
