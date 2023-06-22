library(data.table)
library(readr)

expr <- fread("/scratch/mjpete11/linear_models/data/convert_format/combat_seq_filtered.csv")
expr$V1 <- NULL
expr[1:5,1:5]

small <- expr[1:10,1:10]
small$V1 <- NULL
small[1:5,1:5]
#write_tsv(small, "/scratch/mjpete11/linear_models/data/convert_format/small.tsv")
#write_tsv(expr, "/scratch/mjpete11/linear_models/data/convert_format/combat_seq_filtered.csv")

# Use this file to create the experimental design matrix
organs <- read.csv("/scratch/mjpete11/linear_models/data/organs_combat_seq.csv", sep=",")

# Drop samples from the counts df that are not from paired heart or liver
keeps <- as.vector(organs$SAMPID)
expr1 <- subset(expr, select=keeps) 

# Are there the same number of samples in the expr df and the organs metadata df? 
ncol(expr1)==nrow(organs) # TRUE

# Drop the decimals introduced into the column names
#colnames(expr2) <- keeps

# Reorder the columns (samples) in the expr2 df to be in the
# same order as the samples in the organs metadata df (rows)
# Necessary because voom() assumes the rows of the design matrix are the samples
# and they are in the same order as the column of the expr matrix
idx <- match(organs$SAMPID, colnames(expr1))
expr2 <- as.data.table(expr1)
ordered_expr <- expr2[,..idx]
ordered_expr[1:5,1:5]
expr2[1:5,1:5]
dim(ordered_expr)==dim(expr2) # TRUE TRUE
all(organs$SAMPID==colnames(ordered_expr)) # TRUE
class(ordered_expr) # data.table data.frame
ordered_expr <- as.data.frame(ordered_expr)
ordered_expr1 <- ordered_expr[!duplicated(as.list(ordered_expr))]
ordered_expr1[1:5,1:5]
dim(ordered_expr1)
dim(organs)

# Drop columns in the expression matrix that are not in the design matrix (keep paired samples only)
# Drop samples from the counts df that are not from paired heart or liver
keeps <- as.vector(organs$SAMPID)
ordered_expr2 <- subset(ordered_expr1, select=keeps) 
ncol(ordered_expr2)==nrow(organs) # TRUE

# Drop all columns except the organ and sample columns to use with metabolizer
organs$SUBJID <- NULL
organs$gene <- NULL
organs$combat_seq_counts <- NULL
organs$SMRIN <- NULL
organs$SMTSISCH <- NULL
organs$SMGEBTCH <- NULL
organs$SMGEBTCHT <- NULL

head(organs)
tail(organs)

# Drop duplicate rows
organs_unique <- organs[!duplicated(organs),]
head(organs_unique)
tail(organs_unique)
# Drop samples from the counts df that are not from paired heart or liver
keeps <- as.vector(organs_unique$SAMPID)
ordered_expr3 <- subset(ordered_expr2, select=keeps) 
ncol(ordered_expr3)==nrow(organs_unique) # TRUE
dim(organs_unique) # 296 2
dim(ordered_expr3) # 51259 296
all(organs_unique$SAMPID==colnames(ordered_expr3)) # TRUE

# Write design matrix with only the organ and sample ID
write_tsv(organs_unique, "/scratch/mjpete11/linear_models/data/convert_format/ordered_organs_unique.tsv")

# Expression matrix with the samples in the same order as the design matrix
write_tsv(ordered_expr3, "/scratch/mjpete11/linear_models/data/convert_format/ordered_combat_seq_filtered.tsv")

