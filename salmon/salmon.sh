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

fq1=/scratch/mjpete11/subread_simulations/scripts/single_transcript/100_reads_SLC25A4_R1.fastq
fq2=/scratch/mjpete11/subread_simulations/scripts/single_transcript/100_reads_SLC25A4_R2.fastq
OUTPUT=mapped_reads/SLC25A4.fasta

salmon quant -i $SALMON_INDEX -p $THREADS -l $LIBTYPE -1 $fq1 -2 $fq2 -o $OUTPUT
