#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

BAMS=($(readlink -f slc25_reads/*bam))
item=${BAMS[0]}

# FASTQC on just the slc25 reads
for item in ${BAMS[@]} 
do
	fastqc -o fastqc_results $item
done

# MULTIQC
multiqc fastqc_results

