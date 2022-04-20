#!/bin/bash 
#SBATCH -o sbatch.out
#SBATCH -p htc 

module purge
source /home/mjpete11/.bashrc
conda activate data_science 

# Path
GTEX_PATH=/data/CEM/shared/controlled_access/GTEX/RNA/liver/bam
FILEEND=Aligned.sortedByCoord.out.patched.md

# Define a list
BAMS=($(readlink -f $GTEX_PATH/*bam))
item=${BAMS[0]}

for item in ${BAMS[@]}
do
    outname=$(basename $item)
    outname=slc25_reads/${outname%.bam}_mapped.bam
    samtools view -h -L slc25_genes.bed $item > $outname
done

# FASTQC + MUTLIQC on the BAM files subset to just reads that map to slc25 genes
SUBSET_BAMS=($(readlink -f /scratch/mjpete11/GTEx/liver/slc25_multiqc/slc25_reads/*bam))

for file in ${SUBSET_BAMS[@]}
do
	fastqc -o fastqc_results $file 
done

multiqc .
