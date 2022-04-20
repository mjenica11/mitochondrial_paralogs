#!/bin/bash 
#SBATCH -o sbatch.out
#SBATCH -p htc 

module purge
source /home/mjpete11/.bashrc
conda activate data_science 

# Path
SAM_PATH=/scratch/mjpete11/sweet_data/slc25_multiqc/subset_sams/

# Define a list
SAMS=($(readlink -f $SAM_PATH/*sam))
item=${SAMS[0]}

for item in ${SAMS[@]}
do
    outname=$(basename $item)
    outname=slc25_reads/${outname%.sam}_sapped.bam
    samtools view -h -L slc25_genes.bed $item > $outname
done

# FASTQC + MUTLIQC on the BAM files subset to just reads that map to slc25 genes
#SUBSET_SAMS=($(readlink -f /scratch/mjpete11/GTEx/liver/slc25_multiqc/slc25_reads/*bam))

#for file in ${SUBSET_SAMS[@]}
#do
#	fastqc -o fastqc_results $file 
#done

#multiqc .
