#!/bin/bash

#SBATCH -n 8
#SBATCH -o salmon_index.out
#SBATCH -p general
#SBATCH -t 3-0

module purge
source /home/mjpete11/.bashrc
module load mamba/latest
source activate salmon

#salmon index -t GRCh38_latest_rna.fna -i GRCh38_index
salmon index -t /scratch/mjpete11/mitochondrial_paralogs/sweet_data/salmon_bams/salmon_quantified_sweet_data/index/ncbi_dataset/data/GCF_000001405.40/rna.fna -i GRCh38_index

