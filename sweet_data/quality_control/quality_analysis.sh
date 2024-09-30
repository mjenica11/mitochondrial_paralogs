#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_qc.out
#SBATCH -p general
#SBATCH -t 3-0

module purge
source /home/mjpete11/.bashrc
module mamba/latest
source activate fastqc_env


for i in {784..847}
do
	count=$(printf $i)
	fastqc -o fastqc_results /scratch/mjpete11/mitochondrial_paralogs/sweet_data/salmon_bams/salmon_quantified_sweet_data/salmon_results/SRR7426${count}.fastq
done

multiqc .
