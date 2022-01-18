#!/bin/bash

Rscript -<<EOF
library(stringr)
df<-do.call(cbind,lapply(list.files(path="mapped_reads/",pattern="sf",recursive=T,full.names=T),function(i){
cat(paste(i,"\n"))
read.table(i,header=T)[["NumReads"]]
}))
colnames(df)<-str_extract(list.files(path="mapped_reads/",pattern="sf",recursive=T,full.names=T),pattern="sample_[0-9]+")
write.csv(df,"num_reads.txt",row.names=F)
EOF
