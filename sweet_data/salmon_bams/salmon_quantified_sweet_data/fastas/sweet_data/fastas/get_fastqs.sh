#!/bin/bash
#SBATCH -n 8
#SBATCH -o sbatch_fasta.out
#SBATCH -p general
#SBATCH -t 3-0

module purge
source /home/mjpete11/.bashrc
module load mamba/latest
source activate sra_tools

# Download the FASTQ files from Sweet et al 2018; GSE116250
#fasterq-dump SRR7426784
#fasterq-dump SRR7426785
#fasterq-dump SRR7426786
#fasterq-dump SRR7426787
#fasterq-dump SRR7426788
#fasterq-dump SRR7426789
#fasterq-dump SRR7426790
#fasterq-dump SRR7426791
#fasterq-dump SRR7426792
#fasterq-dump SRR7426793
#fasterq-dump SRR7426794
#fasterq-dump SRR7426795
#fasterq-dump SRR7426796
#fasterq-dump SRR7426797
#fasterq-dump SRR7426798
#fasterq-dump SRR7426799
#fasterq-dump SRR7426800
#fasterq-dump SRR7426801
#fasterq-dump SRR7426802
#fasterq-dump SRR7426803
#fasterq-dump SRR7426804
#fasterq-dump SRR7426805
#fasterq-dump SRR7426806
#fasterq-dump SRR7426807
#fasterq-dump SRR7426808
#fasterq-dump SRR7426809
#fasterq-dump SRR7426810
#fasterq-dump SRR7426811
#fasterq-dump SRR7426812
#fasterq-dump SRR7426813
#fasterq-dump SRR7426814
#fasterq-dump SRR7426815
#fasterq-dump SRR7426816
#fasterq-dump SRR7426817
#fasterq-dump SRR7426818
#fasterq-dump SRR7426819
#fasterq-dump SRR7426820
#fasterq-dump SRR7426821
fasterq-dump SRR7426822
#fasterq-dump SRR7426823
#fasterq-dump SRR7426824
#fasterq-dump SRR7426825
#fasterq-dump SRR7426826
#fasterq-dump SRR7426827
#fasterq-dump SRR7426828
#fasterq-dump SRR7426829
#fasterq-dump SRR7426830
#fasterq-dump SRR7426831
#fasterq-dump SRR7426832
#fasterq-dump SRR7426833
#fasterq-dump SRR7426834
#fasterq-dump SRR7426835
#fasterq-dump SRR7426836
#fasterq-dump SRR7426837
#fasterq-dump SRR7426838
#fasterq-dump SRR7426839
#fasterq-dump SRR7426840
#fasterq-dump SRR7426841
#fasterq-dump SRR7426842
#fasterq-dump SRR7426843
#fasterq-dump SRR7426844
#fasterq-dump SRR7426845
#fasterq-dump SRR7426846
#fasterq-dump SRR7426847
