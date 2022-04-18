#!/bin/bash 
#SBATCH -o sbatch.out

module purge
source /home/mjpete11/.bashrc
conda activate multiqc 

# Subset reads that mapped to the SLC25 family and then perform multiQC
declare -a paralog=(NM_032315.3 NM_207348.3 NM_013386.5 NM_014655.4 
NM_003705.5 NM_017875.4 NM_000387.6 NM_173471.4
NM_00110464 NM_031291.4 NM_021833.5 NM_001151.4
NM_138773.4 NM_001349336.2 NM_031947.4 NM_001271641.2
NM_004277.5 NM_018843.4 NM_014251.3 NM_016612.4
NM_030780.5 NM_033412.4 NM_001330988.2 NM_152707.4
NM_031212.4 NM_014342.4 NM_182556.4 NM_003355.3
NM_003356.4 NM_001191061.2 NM_002635.4 NM_014252.4
NM_001010875.4 NM_030631.4 NM_001039355.3 NM_207117.4
NM_003562.5 NM_001320870.2 NM_001143780.3 NM_001126121.2
NM_012140.5 NM_001034172.4 NM_173637.4 NM_024103.3
NM_001321544.2 NM_031481.3 NM_005984.5 NM_006358.4
NM_001636.4 NM_001012755.5 NM_145305.3 NM_001152.5
NM_001282195.2)

# Define a list
for file in {784..847}; do
	for transcript in ${paralog[@]}; do
		cat subset_sams/SRR7426${file}_mapping_info.csv | grep $transcript > \
				mapped_reads/SRR7426${file}_mapped.sam
	done
done

for i in {784..847}; do
		fastqc -o fastqc_results mapped_reads/SRR7426${count}_mapped.sam
done
