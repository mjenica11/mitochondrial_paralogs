#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate data_science

blastn -query SLC25A1.fa -subject unmapped_SLC25A1.fa -out unmapped_results.out
blastn -query SLC25A1.fa -subject mapped_SLC25A1.fa -out mapped_results.out

