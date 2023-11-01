#!/bin/bash

#SBATCH -n 8
#SBATCH -o index.out

module purge
source /home/mjpete11/.bashrc
conda activate salmon_160

salmon index -t gencode.v43.transcripts.fa.gz -i gencode_index
