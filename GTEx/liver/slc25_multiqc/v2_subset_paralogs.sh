#!/bin/bash 
#SBATCH -o sbatch.out
#SBATCH -p htc 

module purge
source /home/mjpete11/.bashrc
conda activate data_science 

# Path
GTEX_PATH=/data/CEM/shared/controlled_access/GTEX/RNA/liver/bam
FILEEND=Aligned.sortedByCoord.out.patched.md
BEDFILE=slc25_genes.bed

# Define a list
BAMS=($(readlink -f $GTEX_PATH/*bam))
item=${BAMS[0]}
SLC25_SAMS=($(readlink -f slc25_reads/*sam))
item2=${SLC25_SAMS[0]}

for item in ${BAMS[@]}
do
    outname=$(basename $item)
    outname=slc25_reads/${outname%.bam}.sam
    samtools view -h $item $BEDFILE > $outname
done

for item2 in ${SLC25_SAMS[@]} 
do
	fastqc -o fastqc_results $item2
done

multiqc fastqc_results
