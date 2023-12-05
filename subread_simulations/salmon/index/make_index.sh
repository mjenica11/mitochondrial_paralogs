#!/bin/bash

#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

salmon index -t GRCh38_latest_rna.fna -i GRCh38_index
