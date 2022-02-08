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
SAM=hisat_alignment.sam
BAM=hisat_alignment.bam
GFF=GRCh38_latest_genomic.gff
GTF=GRCh38_latest_genomic.gtf
NUMREADS=trans_number_reads.count
UNION_NUMREADS=union_trans_number_reads.count
STRICT_NUMREADS=strict_trans_number_reads.count
NONEMPTY_NUMREADS=nonempty_trans_number_reads.count

# Align and sort reads
#hisat2 -q -x $BASENAME -1 $fq1 -2 $fq2 | samtools sort -o $BAM 
#hisat2 -q -x $BASENAME -1 $fq1 -2 $fq2 | samtools sort -o $SAM 

# Make index
#samtools index $BAM 
#samtools index $SAM 

# Convert GFF file to GTF
#gffread $GFF -T -o $GTF

# Count mapped reads per transcript
#htseq-count -r pos -t CDS -f bam --idattr transcript_id $BAM $GTF > $NUMREADS

# Count mapped reads per transcript and report nonuniquely mapped reads with mode union
htseq-count -r pos -t CDS -f bam --idattr transcript_id -m union --nonunique all $BAM $GTF > $UNION_NUMREADS

# Count mapped reads per transcript and report nonuniquely mapped reads with mode intersection-strict 
htseq-count -r pos -t CDS -f bam --idattr transcript_id -m union --nonunique all $BAM $GTF > $STRICT_NUMREADS

# Count mapped reads per transcript and report nonuniquely mapped reads with mode intersection_nonempty
htseq-count -r pos -t CDS -f bam --idattr transcript_id -m union --nonunique all $BAM $GTF > $NONEMPTY_NUMREADS
