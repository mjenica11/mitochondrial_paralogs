# Split the voom-batch corrected GTEx results into two separate dataframes
# by organ for use with MaREA

library(data.table)
library(readr)
library(dplyr)
library(janitor)

# Read in GTEx results after blocking by batch with voom
#expr <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/batch_voom_qnorm_matrix.csv")
#expr$V1 <- NULL
#expr[1:5,1:5]
#tail(colnames(expr))

# Read in GTEx counts
expr <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")
expr[1:5,1:5]

# Rename the columns to be more descrptive
colnames(expr)[1] <- "Ensembl_ID"
colnames(expr)[2] <- "Hugo_ID"
expr[1:5,1:5]

# Use this file to determine which samples derive from which organ
organs <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/organs_combat_seq.csv", sep=",")
organs[1:5,1:5]

# Drop duplicated rows
#organs2 <- organs[!duplicated(organs$SAMPID),]
#organs2[1:5,1:5]
#dim(organs2) # 296 9

# Drop samples from the counts df that are not from paired heart or liver
keeps <- as.vector(organs$SAMPID)
length(keeps) # 15984 
dim(expr) # 56200 17384 
expr1 <- subset(expr, select=keeps) 
all(colnames(expr1)==keeps)==TRUE # TRUE

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

# Append the Hugo ID column back on the ordered paired expression matrix
ordered_expr3$Hugo_ID <- expr$Hugo_ID
# Move the Hugo_ID column to the front
ordered_expr3 <- ordered_expr3 %>% select(Hugo_ID, everything())
ordered_expr3[1:5,1:5]

# MaREA complains there are duplicated rows...investigate
length(unique(ordered_expr3$Hugo_ID)) # 54592
nrow(ordered_expr3)-length(unique(ordered_expr3$Hugo_ID)) # 1608 duplicated rows
dup <- table(ordered_expr3$Hugo_ID)
any(dup) != 1 # FALSE --> There are no duplicates...
ordered_expr4 <- ordered_expr3 %>% distinct(Hugo_ID, .keep_all = FALSE)
dim(ordered_expr4) # 54592 297 -->  no difference in the number of rows
any(duplicated(ordered_expr3$Hugo_ID))==TRUE # TRUE.. so there are duplicates...
any(duplicated(ordered_expr4$Hugo_ID))==TRUE # FALSE.. so there are duplicates but the same number of rows???
nrow(ordered_expr4)==nrow(ordered_expr3) # FALSE
nrow(ordered_expr4) # 54592
nrow(ordered_expr3) # 56200
dim(ordered_expr4)
ordered_expr4[1:5,]

# Write a small subset of the untransformed paired GTEx samples to test MaREA 
#small <- ordered_expr4[1:20,1:20]
#small[1:5,1:5]
#write_tsv(small, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/convert_format/small_untransformed_GTEx_counts.tsv")

# Write the entire dataset to file
#write_tsv(ordered_expr3, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/convert_format/untransformed_GTEx_counts.tsv")

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
head(organs_unique)  
head(colnames(organs_unique))  
organs_unique[1:5,]  
# Make the rownames into a column
SAMPID <- rownames(ordered_expr_trans)
rownames(ordered_expr_trans) <- NULL
ordered_expr_trans <- cbind(SAMPID, ordered_expr_trans)
head(ordered_expr_trans[1]) 
ordered_expr_trans[1:5,1:5]
# Delete the Hugo_ID row so the rows are the same length in the
# ordered_expr_trans and organs_unique dataframes so it can be concatenated withcbind()
ordered_expr_trans <- ordered_expr_trans[-c(1),]
head(colnames(ordered_expr_trans)) 
dim(ordered_expr_trans) # 296 56201
# Combine the sample x gene matrix with the metadata by
combind <- cbind(ordered_expr_trans, organs_unique, by=c("SAMPID"))
combind[1:5,1:5]
dim(combind) # 296 51263
head(combind$organ)
tail(combind$organ)

# Add the HUGO gene names as column names
expr[1:5,1:5]
names <- expr$Hugo_ID
length(colnames(combind)) # 56204 
length(names) # 56200 ; why are there 4 extra genes in the transposed expression matrix?
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
class(t_heart) # matrix, array
# Convert to dataframe so the floats are encoded as numeric and not as strings
t_heart <- as.data.frame(t_heart)
class(t_heart) # data.frame 
# Convert the rownames into a column
t_heart$hugo_id <- rownames(t_heart)
# Move the hugo ID column to the begining
t_heart <- t_heart %>% select(hugo_id, everything())
# Drop the rownames
rownames(t_heart) <- NULL
t_heart[1:5,1:5]

# Repeat data wrangling with the liver data so it is in a format acceptable with MaREA
t_liver <- t(liver_df)
dim(t_liver) # 51261 148
t_liver[1:5,1:5]
# Drop the organ row
t_liver <- t_liver[!(row.names(t_liver)) %in% c("organ"),] 
# Make the sample ID row into the column names
t_liver <- t_liver %>% row_to_names(1)
class(t_liver) # matrix, array
# Convert to dataframe so the floats are encoded as numeric and not as strings
t_liver <- as.data.frame(t_liver)
class(t_liver) # data.frame 
# Convert the rownames into a column
t_liver$hugo_id <- rownames(t_liver)
# Move the hugo ID column to the begining
t_liver <- t_liver %>% select(hugo_id, everything())
# Drop the rownames
rownames(t_liver) <- NULL
t_liver[1:5,1:5]

# Change the sample names to preserve anonymity before sharing
heart_pseudonames <- sprintf("sample_%s",seq(2:ncol(t_heart)))
head(heart_pseudonames)
colnames(t_heart)[2:ncol(t_heart)] <- heart_pseudonames
t_heart[1:5,1:5]

liver_pseudonames <- sprintf("sample_%s",seq(2:ncol(t_liver)))
head(liver_pseudonames)
colnames(t_liver)[2:ncol(t_liver)] <- liver_pseudonames
t_liver[1:5,1:5]

# Write to file
write_tsv(t_heart, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/pseudonames_GTEx_quantile_voom_heart_df.tsv")
write_tsv(t_liver, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/pseudonames_GTEx_quantile_voom_liver_df.tsv")

