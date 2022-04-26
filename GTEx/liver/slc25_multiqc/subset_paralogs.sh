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
SLC25_SAMS=($(readlink -f slc25_reads/*sam))
item2=${SLC25_SAMS[0]}

for item in ${BAMS[@]}
do
    outname=$(basename $item)
    outname=slc25_reads/${outname%.bam}.sam
    outname=slc25_reads/${outname%.bam}_mapped.bam
    samtools view -h $item chr1:9539465-9585173 chr1:15736258-15741392 \
			chr1:108134043-108200343 chr1:156194104-156212796 \
			chr2:171783405-171894244 chr3:39383370-39397351 \
			chr3:48856926-48898882 chr3:66133610-66378923 \
			chr3:140941836-140980995 chr4:127730400-127774292 \
			chr4:140559431-140568961 chr4:185143266-185150382 \
			chr5:110739007-110765157 chr5:135579172-135888637 \
			chr5:141302635-141304049 chr6:36968141-36986196 \
			chr6:46652975-46678190 chr7:87833568-87876360 \
			chr7:96120220-96322098 chr8:23528956-23575463 \
			chr8:103398638-103415107 chr9:37887604-37904127 \
			chr9:128068232-128109245 chr10:68482344-68527523 \
			chr10:99610522-99620439 chr11:790475-797998 \
			chr11:47617317-47642559 chr11:65375192-65382646 \
			chr11:73974672-73982843 chr11:74000277-74009085 \
			chr12:98593686-98606367 chr13:40789611-40812460 \
			chr13:45393316-45418373 chr14:36677921-37172606 \
			chr14:100291116-100306444 chr14:100323339-100330421 \
			chr17:4937130-4940046 chr17:8290061-8295400 \
			chr17:44319628-44324823 chr17:75272992-75289433 \
			chr17:81712284-81721012 chr18:31759584-31760880 \
			chr19:6426037-6433763 chr19:6440064-6459792 \
			chr19:19063994-19113030 chr22:17563498-17590995 \
			chr22:40769630-40819346 chr22:219175581-19178736 \
			chrX:1386152-1392113 chrX:104099214-104157009 \
			chrX:119399336-119454478 chrX:119468444-119471396 \
			chrX:130339919-130373357 > $outname
    samtools view -h -b $item slc25_genes.bed > $outname
	samtools index $item2
done

for item2 in ${SLC25_SAMS[@]} 
do
	fastqc -o fastqc_results $item2
done

multiqc fastqc_results
