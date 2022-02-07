#!/bin/bash

#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

# Manually downloaded RefSeq transcripts from:
# https://www.ncbi.nlm.nih.gov/genome/guide/human/

# Build HGFM index
hisat2-build -p 16 GRCh38_latest_rna.fna GRCh38
