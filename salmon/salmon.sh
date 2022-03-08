#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate salmon_160

SALMON_INDEX=/scratch/mjpete11/salmon/index/GRCh38_index/
LIBTYPE=A
THREADS=8

for i in {1..53}
do
		count=$(printf "%02d" $i)
		fq1=/scratch/mjpete11/subread_simulations/scripts/individual_simulations/SLC25A${count}_R1.fastq
		fq2=/scratch/mjpete11/subread_simulations/scripts/individual_simulations/SLC25A${count}_R2.fastq
		OUTPUT=individual_simulations/SLC25A${count}
		MAPPING=individual_simulations/SLC25A${count}/mapping_info.csv

		salmon quant -i $SALMON_INDEX -p $THREADS -l $LIBTYPE -1 $fq1 -2 $fq2 \
			--dumpEq --validateMappings --writeMappings $MAPPING \
			--writeUnmappedNames -o $OUTPUT
done
