library(stringr)
library(Biostrings)

# Read in the fasta the reads were simulated from
fa<-readAAStringSet("/scratch/mjpete11/subread_simulations/fastas/SLC25A4.fa")

# Get the ENST* IDs from the fasta the reads were simulated from
ENST_names<-sapply(strsplit(names(fa),"\\|"),"[[",1)

# Drop the extra space at the end of the IDs
ENST_names <- trimws(ENST_names,which=c("right"),whitespace="[ \t\r\n]")
ENST_names

# Get the ENST* IDs from the fasta the reads were simulated from
SLC25_names<-sapply(strsplit(names(fa),"\\|"),"[[",2)

# Drop the extra space at the end of the IDs
SLC25_names <- trimws(SLC25_names,which=c("left"),whitespace="[ \t\r\n]")
SLC25_names

# Extract the NumReads column from the salmon quant.sf files
df<-do.call(cbind,
            lapply(
			list.files(path="test_mapped_reads/",
					   pattern="sf",
					   recursive=T,
					   full.names=T),
			function(i){
				cat(paste(i,"\n"))
				read.table(i,header=T)[["NumReads"]]
			}))

# Add ENST* IDs as column names and rownames
#colnames(df)<-ENST_names
colnames(df)<-"NumReads"

# Extract the transcript IDs from one quant.sf file and append to the df
df2 <- read.csv("test_mapped_reads/100_reads_SLC25A4.fasta/quant.sf", sep="\t")
transcript_IDs <- df2$Name
df <- as.data.frame(cbind(transcript_IDs, df))

# Write the same df except use the SLC25* nomenclature for the colnames
df3 <- df
#colnames(df3)[2:ncol(df3)] <- SLC25_names
colnames(df3)[2:ncol(df3)] <- "sample1" 

# Write to csv
write.csv(df,"output/RefSeq_num_reads.txt",row.names=F)
#write.csv(df3,"output/SLC25_num_reads.txt",row.names=F)



#Alternative way to write the column names with the SLC25* IDs 
#colnames(df)<-str_extract(list.files(path="mapped_reads/",pattern="sf",recursive=T,full.names=T),pattern="SLC25A[0-9]+")
# Code to scale up to n samples
#colnames(df)<-str_extract(list.files(path="mapped_reads/",pattern="sf",recursive=T,full.names=T),pattern="sample_[0-9]+")
