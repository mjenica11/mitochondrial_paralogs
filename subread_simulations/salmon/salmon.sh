#!/bin/bash
#SBATCH -o sbatch_test.out
#SBATCH -t 0-00:48:00   # time in d-hh:mm:ss

module purge
module load salmon-1.4.0-gcc-11.2.0

SALMON_INDEX=/scratch/mjpete11/mitochondrial_paralogs/salmon/index/GRCh38_index/
LIBTYPE=A
THREADS=8

for i in {1..9}
do
		count=$(printf "%02d" $i)
		fq1=/scratch/mjpete11/mitochondrial_paralogs/simulations/1000_simulations/sample_0${count}_1.fasta
		fq2=/scratch/mjpete11/mitochondrial_paralogs/simulations/1000_simulations/sample_0${count}_2.fasta
		OUTPUT=error_5_individual_simulations/sample_${count}
		MAPPING=error_5_individual_simulations/sample_${count}/mapping_info.csv

		salmon quant -i $SALMON_INDEX -p $THREADS -l $LIBTYPE -1 $fq1 -2 $fq2 \
			--dumpEq --validateMappings --writeMappings $MAPPING \
			--writeUnmappedNames -o $OUTPUT
done
