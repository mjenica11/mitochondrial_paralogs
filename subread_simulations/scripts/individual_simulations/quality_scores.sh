#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate fastqc_environment

for i in {1..53}
do
	count=$(printf "%02d" $i)
	FQ1=SLC25A${count}_R1.fastq
	FQ2=SLC25A${count}_R2.fastq
	fastqc -o fastqc_results $FQ1 $FQ2 
done

multiqc .
