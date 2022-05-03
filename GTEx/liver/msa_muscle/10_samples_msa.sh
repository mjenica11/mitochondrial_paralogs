#!/bin/bash
#SBATCH -o sbatch.out

module purge
source /home/mjpete11/.bashrc
conda activate clustalo

FILEPATH=/scratch/mjpete11/GTEx/liver/msa_muscle

# Align 10 reads that mapped to SLC25A18 and the mRNA sequence they were
# simulated from
for file in {1..10}
do
	muscle -in $FILEPATH/mapped_fastas/mapped_GTEX-ZZPU-0426-SM-5GZYH_${file}.fa \
			-clwstrictout $FILEPATH/results/mapped_GTEX-ZZPU-0426-SM-5GZYH_${file}.clw
done
