#!/bin/bash 
#SBATCH -o sbatch.out
#SBATCH -p phi

module purge
source /home/mjpete11/.bashrc
conda activate samtools 

# List of GTEx sample IDs
array=(11DXY-0526-SM-5EGGQ 11EQ9-0526-SM-5A5JZ 11DXZ-0126-SM-5EGGY)

# Path
PATH=/data/CEM/shared/controlled_access/GTEX/RNA/liver/bam
FILEEND=Aligned.sortedByCoord.out.patched.md

# Define a list
for item in ${!array[@]}; do
	samtools view -h -b -L slc25_genes.bed $PATH/GTEX-${array[$item]}.$FILEEND.bam /
	>> slc25_reads/GTEX-${array[$item]}_mapped.bam
done

#for i in ${!array}; do
#		fastqc -o fastqc_results mapped_reads/GTEX-${array[$index]}_mapped.bam
#done
