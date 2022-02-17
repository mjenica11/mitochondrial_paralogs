#!/bin/bash
#SBATCH -n 8
#SBATCH -o hisat_test.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate hisat2 

fq1=/scratch/mjpete11/subread_simulations/scripts/output/1000_reads_SLC25A4_R1.fastq
fq2=/scratch/mjpete11/subread_simulations/scripts/output/1000_reads_SLC25A4_R2.fastq
BASENAME=GRCh38
SAM=hisat_alignment.sam
BAM=hisat_alignment.bam
GFF=GRCh38_latest_genomic.gff
GTF=GRCh38_latest_genomic.gtf

# Align and sort reads
# Generate a BAM file to use with IGV viewer and human readable SAM file
hisat2 -q -x $BASENAME -1 $fq1 -2 $fq2 | samtools sort -o $BAM 
hisat2 -q -x $BASENAME -1 $fq1 -2 $fq2 | samtools sort -o $SAM 

# Make index
# Generate a BAM index to use with IGV viewer and a human readable SAM file
samtools index $BAM 
samtools index $SAM 

# Convert GFF file to GTF
#gffread $GFF -T -o $GTF
