#!/bin/bash
# Convert sam files of SLC25 reads to counts for mediation analysis

# Path to heart SLC counts: /scratch/mjpete11/GTEx/heart/slc25_reads
# Path to liver SLC counts: /scratch/mjpete11/GTEx/liver/slc25_multiqc/slc25_reads
HEART=/scratch/mjpete11/GTEx/heart/slc25_reads

cut -f 10 $HEART/GTEX-18A66-1226-SM-7LT8G.Aligned.sortedByCoord.out.patched.md.sam \
		| sort | uniq -c > heart.txt 
