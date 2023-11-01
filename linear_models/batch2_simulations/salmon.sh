#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out

module purge
source /home/mjpete11/.bashrc
conda activate salmon_160

SALMON_INDEX=/scratch/mjpete11/linear_models/batch2_simulations/gencode_index/
LIBTYPE=A
THREADS=8

for i in {299..405}
do
		#count=$(printf "%02d" $i)
		count=$i
		fq1=/scratch/mjpete11/linear_models/batch2_simulations/1000_samples/sample_${count}_1.fasta
		fq2=/scratch/mjpete11/linear_models/batch2_simulations/1000_samples/sample_${count}_2.fasta
		OUTPUT=individual_simulations/sample_${count}
		MAPPING=individual_simulations/sample_${count}/mapping_info.csv

		salmon quant -i $SALMON_INDEX -p $THREADS -l $LIBTYPE -1 $fq1 -2 $fq2 \
			--minAssignedFrags 1 --dumpEq --validateMappings --writeMappings $MAPPING \
			--writeUnmappedNames -o $OUTPUT
done
