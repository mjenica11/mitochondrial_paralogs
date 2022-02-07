#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate hisat2 

fq1=/scratch/mjpete11/subread_simulations/scripts/single_transcript/100_reads_SLC25A4_R1.fastq
fq2=/scratch/mjpete11/subread_simulations/scripts/single_transcript/100_reads_SLC25A4_R2.fastq
#HISAT2_INDEXES=/scratch/mjpete11/hisat/index
BASENAME=GRCh38
SAM=/scratch/mjpete11/hisat/hisat_alignment.sam

hisat2 -q -x $BASENAME -1 $fq1 -2 $fq2 | samtools view -o $SAM
