#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate fastqc_environment

for i in {784..847}
do
	count=$(printf $i)
	fastqc -o fastqc_results SRR7426${count}.fastq
done

multiqc .
