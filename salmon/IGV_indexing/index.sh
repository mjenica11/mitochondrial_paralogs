#!/bin/bash
#SBATCH -o sbatch.out

# Index the salmon qualitification file of the simulated SLC25A51/52 reads
# in order to view with IGView ;
# Copied the mapping_info.csv file from salmon/individual_simulations/SLC25A52
# and renames it to SLC25A52.sam and repeated for SLC25A51
module purge
source /home/mjpete11/.bashrc
conda activate data_science 

# Convert sam to bam format; recommended format for IGViewer
samtools view -h -o SLC25A51.bam SLC25A51.sam 
samtools view -h -o SLC25A52.bam SLC25A52.sam 

# Sort and index
samtools sort SLC25A52.bam -o SLC25A52.bam.bai
samtools index SLC25A52.bam.bai prepared_SLC25A52.bam.bai

samtools sort SLC25A51.bam -o SLC25A51.bam.bai
samtools index SLC25A51.bam.bai prepared_SLC25A51.bam.bai
