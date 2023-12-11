# Split the voom-batch corrected GTEx results into two separate dataframes
# by organ for use with MaREA

library(data.table)
library(readr)
library(dplyr)
library(janitor)

# Read in GTEx results after blocking by batch with voom
expr <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/batch_voom_qnorm_matrix.csv")
expr$V1 <- NULL
expr[1:5,1:5]
tail(colnames(expr))

# Use this file to determine which samples derive from which organ
organs <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/organs_combat_seq.csv", sep=",")

# Drop samples from the counts df that are not from paired heart or liver
keeps <- as.vector(organs$SAMPID)
expr1 <- subset(expr, select=keeps) 

# Are there the same number of samples in the expr df and the organs metadata df? 
ncol(expr1)==nrow(organs) # TRUE

# Drop the decimals introduced into the column names
#colnames(expr2) <- keeps

# Reorder the columns (samples) in the expr2 df to be in the
# same order as the samples in the organs metadata df (rows)
# This should make it easier to split into two dataframes by organ type
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

# Drop all columns except the organ and sample columns to use with MaREA 
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
ordered_expr3[1:5,1:5] 
dim(ordered_expr3) # 51259 296
all(organs_unique$SAMPID==colnames(ordered_expr3)) # TRUE
head(organs_unique) # 296 2
tail(organs_unique) # 296 2

# Write design matrix with only the organ and sample ID
#write_tsv(organs_unique, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/convert_format/ordered_organs_unique.tsv")

# Transpose the expression matrix into a sample x gene format
ordered_expr_trans <- t(ordered_expr3)
ordered_expr_trans[1:5,1:5]
# Append the organ column from the organs_unique df to the expression matrix
# by the sample ID
dim(ordered_expr_trans) # 296 51259
class(ordered_expr_trans) # matrix
ordered_expr_trans <- as.data.frame(ordered_expr_trans)
class(ordered_expr_trans) # data.frame 
ordered_expr_trans[1:3,] 
head(rownames(ordered_expr_trans)) 
head(colnames(ordered_expr_trans)) 
dim(organs_unique) # 296 2 
head(organs_unique) # 296 2 
organs_unique[1:5,]  
# Make the rownames into a column
SAMPID <- rownames(ordered_expr_trans)
rownames(ordered_expr_trans) <- NULL
ordered_expr_trans <- cbind(SAMPID, ordered_expr_trans)
head(ordered_expr_trans[1]) 
head(ordered_expr_trans) 
# Combine the sample x gene matrix with the metadata by
combind <- cbind(ordered_expr_trans, organs_unique, by=c("SAMPID"))
combind[1:5,1:5]
dim(combind) # 296 51263
head(combind$organ)
tail(combind$organ)

# Add the HUGO gene names as column names
expr[1:5,1:5]
names <- expr$hugo_names
length(colnames(combind)) # 51263
length(names) # 51259 ; why are there 4 extra genes in the transposed expression matrix?
####### there are 2 SAMPID column, organ col, and a "by" column...
combind[1] <- NULL
head(combind$by)
combind$by <- NULL
# Move the SAMPID ID and organ column to the front
combind <- combind %>% select(SAMPID, everything())
combind[1:5,1:5]
combind <- combind %>% select(organ, everything())
combind[1:5,1:5]
colnames(combind[2:ncol(combind)]) <- names
combind[1:5,1:5]
head(colnames(combind))
dim(combind) # 296 51261
# Append the hugo gene names as the column names
length(names)==length(colnames(combind))-2 # TRUE
head(colnames(combind)[3:ncol(combind)])
colnames(combind)[3:ncol(combind)] <- names
head(colnames(combind)[3:ncol(combind)])
combind[1:5,1:5]

# Split into two separate dataframes
split_df <- split(combind, combind$organ, drop=FALSE)
class(split_df) # list
str(split_df) # heart_df and liver_df
split_df$liver[1:5,1:5]
split_df$heart[1:5,1:5]
# Store resulting dataframes and separate variables
heart_df <- split_df$heart
head(heart_df)
class(heart_df)
liver_df <- split_df$liver
head(liver_df)
class(liver_df)
# Transpose heart and liver dataframes so they are in a gene x sample format
# for use with MaREA
t_heart <- t(heart_df)
dim(t_heart) # 51261 148
t_heart[1:5,1:5]
# Drop the organ row
t_heart <- t_heart[!(row.names(t_heart)) %in% c("organ"),] 
# Make the sample ID row into the column names
t_heart <- t_heart %>% row_to_names(1)
# Convert the values from character to float
t_heart <- sapply(t_heart, as.numeric)
###### stopped here
###### idk why the next step causes it to fail
# Convert the rownames into a column
t_heart$hugo_id <- rownames(t_heart)
t_heart[1:5,1:5]
t_liver <- t(liver_df)
dim(t_liver) # 51261 148


write_tsv(heart_df, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/convert_format/heart_df.tsv")
write_tsv(liver_df, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/convert_format/liver_df.tsv")

