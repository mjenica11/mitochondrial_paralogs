#!/bin/bash
#SBATCH -o sbatch.out
#SBATCH -p phi

module purge
source /home/mjpete11/.bashrc
conda activate data_science

# Subset reads that mapped to SLC25A1
samtools view -h -b -L SLC25A1.bed GTEX-ZZPU-0426-SM-5GZYH.Aligned.sortedByCoord.out.patched.md.bam \
		> SLC25A1_GTEX-ZZPU-0426-SM-5GZYH.bam

# Convert bam to sam
samtools view -h -o SLC25A1_GTEX-ZZPU-0426-SM-5GZYH.sam SLC25A1_GTEX-ZZPU-0426-SM-5GZYH.bam
