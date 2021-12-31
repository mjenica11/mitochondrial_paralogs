#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out

conda activate salmon_160
SALMON_INDEX=/scratch/mjpete11/salmon_alignment/default
LIBTYPE=A
threads=8

for i in {1..20}
do
    count=$(printf "%02d" $i)
    fq1=simulated_reads/sample_${count}_1.fasta
    fq2=simulated_reads/sample_${count}_2.fasta
    OUTPUT=mapped_reads/sample_${count}

    
    salmon quant -i $SALMON_INDEX -p $threads -l $LIBTYPE -1 $fq1 -2 $fq2 -o $OUTPUT
done
