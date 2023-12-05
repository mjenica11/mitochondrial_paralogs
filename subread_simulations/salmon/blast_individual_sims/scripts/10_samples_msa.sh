#!/bin/bash
#SBATCH -o sbatch.out

module purge
source /home/mjpete11/.bashrc
conda activate clustalo

FILE_PATH=/scratch/mjpete11/salmon/blast_individual_sims

# align mapped reads of a non-failing heart control with muscle
for file in {1..10}
do
	muscle -in $FILE_PATH/mapped_fastas/mapped_SLC25A1_${file}.fa \
			-clwstrictout  $FILE_PATH/mapped_results/SLC25A1_${file}.clw
done
