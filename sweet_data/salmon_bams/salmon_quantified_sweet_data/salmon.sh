#!/bin/bash
#SBATCH -n 8
#SBATCH -o bams_attempt5.out
#SBATCH -p general
#SBATCH -t 3-0

module purge
source /home/mjpete11/.bashrc
module load mamba/latest
source activate salmon

SALMON_INDEX=/scratch/mjpete11/mitochondrial_paralogs/sweet_data/salmon_bams/salmon_quantified_sweet_data/index/GRCh38_index/
LIBTYPE=A
THREADS=8

for i in {784..847}
do
	count=$(printf $i)
	FASTQ=/scratch/mjpete11/mitochondrial_paralogs/sweet_data/salmon_bams/salmon_quantified_sweet_data/fastas/sweet_data/fastas/SRR7426${count}.fastq
	FASTQ=/scratch/mjpete11/mitochondrial_paralogs/sweet_data/salmon_bams/salmon_quantified_sweet_data/fastas/sweet_data/fastas/SRR7426${count}.fastq
	OUTPUT=./salmon_results/SRR7426${count}
	MAPPING=./salmon_results/SRR7426${count}/mapping_info.csv

	salmon quant -i $SALMON_INDEX -p $THREADS -l $LIBTYPE -r $FASTQ \
			     --dumpEq --validateMappings --writeMappings $MAPPING \
				 --writeUnmappedNames -o $OUTPUT
#done
