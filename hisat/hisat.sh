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
GFF=GRCh38_latest_genomic.gff
GTF=GRCh38_latest_genomic.gtf
NUMREADS=number_reads.count

# Align and sort reads
hisat2 -q -x $BASENAME -1 $fq1 -2 $fq2 | samtools sort -o $BAM 

# Make index
samtools index $BAM 

# Convert GFF file to GTF
gffread $GFF -T -o $GTF

# Count mapped reads per transcript
htseq-count -r pos -t CDS -f bam $BAM $GTF > $NUMREADS
