#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate salmon_160

SALMON_INDEX=/scratch/mjpete11/salmon/default
LIBTYPE=A
THREADS=8

fq1=/scratch/mjpete11/subread_simulations/scripts/output/combined_R1.fastq
fq2=/scratch/mjpete11/subread_simulations/scripts/output/combined_R2.fastq
OUTPUT=mapped_reads/combined.fasta

salmon quant -i $SALMON_INDEX -p $THREADS -l $LIBTYPE -1 $fq1 -2 $fq2 -o $OUTPUT

# Code to quantify n samples
#for i in {1..1000}
#do
#   count=$(printf "%02d" $i)
#	fq1=/scratch/mjpete11/subread_simulations/scripts/sample_${count}_R1.fasta
#	fq1=/scratch/mjpete11/subread_simulations/scripts/sample_${count}_R2.fasta
#	OUTPUT=mapped_reads/$sample_{count}.fasta 
    
#	salmon quant -i $SALMON_INDEX -p $threads -l $LIBTYPE -1 $fq1 -2 $fq2 -o $OUTPUT
#done
