#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate data_science

# Count the number of unmapped reads in the GTEx liver samples 
#samtools view -c -f 4 /data/CEM/shared/controlled_access/GTEX/RNA/liver/bam/GTEX-ZF2S-3026-SM-4WWCH.Aligned.sortedByCoord.out.patched.md.bam \
#		| cut -f1 > GTEX-ZF2S-3026-SM-4WWCH_unmapped_reads.txt

# Count the total number of reads in the GTEx liver samples
#samtools view -c /data/CEM/shared/controlled_access/GTEX/RNA/liver/bam/GTEX-ZF2S-3026-SM-4WWCH.Aligned.sortedByCoord.out.patched.md.bam \
#		| cut -f1 > GTEX-ZF2S-3026-SM-4WWCH_total_reads.txt

# Divide the length of the total number of reads by the unmapped reads
nom=$(cat GTEX-ZF2S-3026-SM-4WWCH_unmapped_reads.txt) 
denom=$(cat GTEX-ZF2S-3026-SM-4WWCH_total_reads.txt)
res=$((nom/denom))
echo "$res" >> liver_read_proportions.txt
