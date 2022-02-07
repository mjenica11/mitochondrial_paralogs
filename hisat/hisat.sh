#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate hisat2 

fq1=/scratch/mjpete11/subread_simulations/scripts/single_transcript/100_reads_SLC25A4_R1.fastq
fq2=/scratch/mjpete11/subread_simulations/scripts/single_transcript/100_reads_SLC25A4_R2.fastq
BASENAME=GRCh38
BAM=hisat_alignment.bam

hisat2 -q -x $BASENAME -1 $fq1 -2 $fq2 | samtools sort -o $BAM 
samtools index $BAM 
