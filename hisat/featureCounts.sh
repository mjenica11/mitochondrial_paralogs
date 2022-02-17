#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_test.out
#SBATCH -p htc

module purge
source /home/mjpete11/.bashrc
conda activate hisat2 

BAM=hisat_alignment.bam
GTF=GRCh38_latest_genomic.gtf
NEW_GTF=corrected.gtf
NUMREADS=featureCounts.count

# featureCounts want the gene_ID column to be the 9th column
# Convert gene_ID column to be the 9th column in the GTF file
perl -ne 'chomp; @a=split/\t/; %h=split(/ /,$a[8]); 
$a[8]=join(" ",("gene_id",$h{"gene_id"},"transcript_id",$h{"transcript_id"}));
print join("\t",@a),"\n";' $GTF > corrected.gtf

# Count the number of reads mapped per feature with featureCounts
featureCounts -a $NEW_GTF -o $NUMREADS $BAM
