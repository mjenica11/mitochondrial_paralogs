#!/bin/bash
#SBATCH -n 8
#SBATCH -o bams.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate salmon_160

SALMON_INDEX=/scratch/mjpete11/salmon/index/GRCh38_index/
LIBTYPE=A
THREADS=8

for i in {784..847}
do
		count=$(printf $i)
		FASTQ=/scratch/mjpete11/sweet_data/SRR7426${count}.fastq
		OUTPUT=SRR7426${count}
		MAPPING=SRR7426${count}/mapping_info.csv

		salmon quant -i $SALMON_INDEX -p $THREADS -l $LIBTYPE -r $FASTQ \
			--dumpEq --validateMappings --writeMappings $MAPPING \
			--writeUnmappedNames -o $OUTPUT
done
