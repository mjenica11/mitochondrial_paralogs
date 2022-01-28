#!/bin/bash

#SBATCH -q normal 
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=%u@asu.edu

module purge
source /home/mjpete11/.bashrc
conda activate data_science

esearch -db nucleotide -query "NM_032315.3" | efetch -format fasta > SLC25A33.fa
esearch -db nucleotide -query "NM_207348.3" | efetch -format fasta > SLC25A34.fa
esearch -db nucleotide -query "NM_013386.5" | efetch -format fasta > SLC25A24.fa
esearch -db nucleotide -query "NM_014655.4" | efetch -format fasta > SLC25A44.fa
esearch -db nucleotide -query "NM_003705.5" | efetch -format fasta > SLC25A12.fa
esearch -db nucleotide -query "NM_017875.4" | efetch -format fasta > SLC25A38.fa
esearch -db nucleotide -query "NM_000387.6" | efetch -format fasta > SLC25A20.fa
esearch -db nucleotide -query "NM_173471.4" | efetch -format fasta > SLC25A26.fa
esearch -db nucleotide -query "NM_001104647.3" | efetch -format fasta > SLC25A36.fa
esearch -db nucleotide -query "NM_031291.4" | efetch -format fasta > SLC25A31.fa
esearch -db nucleotide -query "NM_021833.5" | efetch -format fasta > SLC25A7.fa
esearch -db nucleotide -query "NM_001151.4" | efetch -format fasta > SLC25A4.fa
esearch -db nucleotide -query "NM_138773.4" | efetch -format fasta > SLC25A46.fa
esearch -db nucleotide -query "NM_001349336.2" | efetch -format fasta > SLC25A48.fa
esearch -db nucleotide -query "NM_031947.4" | efetch -format fasta > SLC25A2.fa
esearch -db nucleotide -query "NM_001271641.2" | efetch -format fasta > SLC25A49.fa
esearch -db nucleotide -query "NM_004277.5" | efetch -format fasta > SLC25A27.fa
esearch -db nucleotide -query "NM_018843.4" | efetch -format fasta > SLC25A40.fa
esearch -db nucleotide -query "NM_014251.3" | efetch -format fasta > SLC25A13.fa
esearch -db nucleotide -query "NM_016612.4" | efetch -format fasta > SLC25A37.fa
esearch -db nucleotide -query "NM_030780.5" | efetch -format fasta > SLC25A32.fa
esearch -db nucleotide -query "NM_033412.4" | efetch -format fasta > SLC25A51.fa
esearch -db nucleotide -query "NM_001330988.2" | efetch -format fasta > SLC25A25.fa
esearch -db nucleotide -query "NM_152707.4" | efetch -format fasta > SLC25A16.fa
esearch -db nucleotide -query "NM_031212.4" | efetch -format fasta > SLC25A28.fa
esearch -db nucleotide -query "NM_014342.4" | efetch -format fasta > SLC25A50.fa
esearch -db nucleotide -query "NM_182556.4" | efetch -format fasta > SLC25A45.fa
esearch -db nucleotide -query "NM_003355.3" | efetch -format fasta > SLC25A8.fa
esearch -db nucleotide -query "NM_003356.4" | efetch -format fasta > SLC25A9.fa
esearch -db nucleotide -query "NM_001191061.2" | efetch -format fasta > SLC25A22.fa
esearch -db nucleotide -query "NM_002635.4" | efetch -format fasta > SLC25A3.fa
esearch -db nucleotide -query "NM_014252.4" | efetch -format fasta > SLC25A15.fa
esearch -db nucleotide -query "NM_001010875.4" | efetch -format fasta > SLC25A30.fa
esearch -db nucleotide -query "NM_030631.4" | efetch -format fasta > SLC25A21.fa
esearch -db nucleotide -query "NM_001039355.3" | efetch -format fasta > SLC25A29.fa
esearch -db nucleotide -query "NM_207117.4" | efetch -format fasta > SLC25A47.fa
esearch -db nucleotide -query "NM_003562.5" | efetch -format fasta > SLC25A11.fa
esearch -db nucleotide -query "NM_001320870.2" | efetch -format fasta > SLC25A35.fa
esearch -db nucleotide -query "NM_001143780.3" | efetch -format fasta > SLC25A39.fa
esearch -db nucleotide -query "NM_001126121.2" | efetch -format fasta > SLC25A19.fa
esearch -db nucleotide -query "NM_012140.5" | efetch -format fasta > SLC25A10.fa
esearch -db nucleotide -query "NM_001034172.4" | efetch -format fasta > SLC25A52.fa
esearch -db nucleotide -query "NM_173637.4" | efetch -format fasta > SLC25A41.fa
esearch -db nucleotide -query "NM_024103.3" | efetch -format fasta > SLC25A23.fa
esearch -db nucleotide -query "NM_001321544.2" | efetch -format fasta > SLC25A42.fa
esearch -db nucleotide -query "NM_031481.3" | efetch -format fasta > SLC25A18.fa
esearch -db nucleotide -query "NM_005984.5" | efetch -format fasta > SLC25A1.fa
esearch -db nucleotide -query "NM_006358.4" | efetch -format fasta > SLC25A17.fa
esearch -db nucleotide -query "NM_001636.4" | efetch -format fasta > SLC25A6.fa
esearch -db nucleotide -query "NM_001012755.5" | efetch -format fasta > SLC25A53.fa
esearch -db nucleotide -query "NM_145305.3" | efetch -format fasta > SLC25A43.fa
esearch -db nucleotide -query "NM_001152.5" | efetch -format fasta > SLC25A5.fa
esearch -db nucleotide -query "NM_001282195.2" | efetch -format fasta > SLC25A14.fa
