#!/bin/bash 
#SBATCH -o sbatch.out
#SBATCH -p htc 

module purge
source /home/mjpete11/.bashrc
conda activate data_science 

# Path
GTEX_PATH=/data/CEM/shared/controlled_access/GTEX/RNA/heart_left_ventricle/bam
FILEEND=Aligned.sortedByCoord.out.patched.md

# Define a list
BAMS=($(readlink -f $GTEX_PATH/*bam))
item=${BAMS[0]}
SLC25_BAMS=($(readlink -f slc25_reads/*bam))
item2=${SLC25_BAMS[0]}

for item in ${BAMS[@]}
do
    outname=$(basename $item)
    outname=slc25_reads/${outname%.bam}_mapped.bam
    samtools view -h -b slc25_genes.bed $item > $outname
	samtools index $item2
	fastqc -o fastqc_results $item2
done

multiqc fastqc_results
