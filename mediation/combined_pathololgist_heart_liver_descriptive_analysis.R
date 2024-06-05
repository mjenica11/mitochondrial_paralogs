# Aim 2 descriptive analysis of combined heart + liver pathOlogist activity scores
# Examine correlation between SLC25 and TCA cycle activity scores

# Library
library(data.table)
library(janitor)
library(tidyverse)
library(dplyr)
library(stringr)
library(mediation)
library(purrr)
library(broom)
library(patchwork) # to make the barplots

# Read in the pathways and activity scores at 50% activity
activity <- fread("~/Downloads/filtered_activity_scores_50percent.csv",sep=",")
activity[1:5,1:5]

any(duplicated(activity))==TRUE # FALSE

# Read in the gene counts with the RIN and activity scores
organs <- read.csv("~/Downloads/qnorm_voom_SLC25_biomarkers_manifest.csv", sep=",")
colnames(organs)[10]<-"SampleNum"
organs[1:5,1:5]
dim(organs) # 34632 10 --> 53 SL25 genes + 64 metabolic biomarkers * 148 samples * 2 samples per individual

# Check for any duplicated rows
any(duplicated(organs))==TRUE # FALSE

# Change the .s to  _ in the sample names in the activity dataframe
colnames(activity) <- gsub(x=colnames(activity), pattern="\\.", replacement="-")
activity[1:5,1:5]

# List of SLC25 paralogs
SLC <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6","UCP1", "UCP2", "UCP3", "SLC25A10", "SLC25A11", "SLC25A12","SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17","SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22","SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27","SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32","SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37","SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42","SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47","SLC25A48", "MTCH1", "MTCH2", "SLC25A51", "SLC25A52", "SLC25A53")

# List of metabolic biomarkers 
biomarkers <- c("SLC16A1", "SLC16A2", "HK1", "HK2", "HK3", "GPI", "PFKM", "PFKP", "ALDOA", "ALDOB", "TPI1", "GAPDH", "PGK1","PGAM1", "PGAM5", "PGAM4", "ENO1","PKLR", "LDHA", "LDHB", "LDHC", "LDHD", "CPT1A","CPT1B", "CPT2", "ACADVL", "ACADL", "ACADM", "ACADS","ACADSB", "ACAD11", "ACAD8", "ACAD10", "ACAD9", "ECHS1","ECH1", "ECI1", "ECI2", "HADHA", "PDHA1", "PDHA2","PDHB", "PDP1", "PDHX", "CS", "ACO1", "ACO2", "IDH3G","IDH1", "IDH2", "IDH3B", "OGDH", "SUCLG2","SUCLG1", "SUCLA2", "SDHB", "SDHD", "SDHA", "SDHC", "SDHAF2","FH", "MDH1", "MDH2", "MTND1")

# Combine list
genes <- c(SLC, biomarkers)
genes
length(genes) # 117

# reshape data from long to wide
activity <- activity %>% reshape2::melt(id.vars=c('pathway_name'))
dim(activity) # 86136 3
activity[1:5,]
head(activity)

# Clearer column names
colnames(activity) <- c("pathway","sample_ID","activity")
class(activity$pathway)
class(activity)

# Subset just the TCA cycle pathway activity scores since it is the dependent variable of interest
TCA_DF <- activity %>% filter(str_detect(pathway,"tca"))
head(TCA_DF)

# What does the TCA cycle activity look like?
plot(density(TCA_DF$activity))

# Stratify the scores into low, medium, and high
quantile(TCA_DF$activity)
summary(TCA_DF$activity)
TCA_DF$Strata <- cut(TCA_DF$activity,breaks = quantile(TCA_DF$activity),labels = paste0("Q",1:4),include.lowest = T)

# Make a dataframe to make a heatmap; add the TCA cycle activity score quartiles
# left_join() keeps all observations in x; keep all observations in the 
# organs dataframe; add the TCA cycle activity score and quartile strata
# by merging on the sample_ID
HeatMAP <- organs %>% left_join(TCA_DF %>% dplyr::select(sample_ID, activity, Strata))
dim(HeatMAP) # 34632 13

# Make sure the number of rows is identical (same number of genes)
# 53 SLC25 + 65 biomarkers * 148 individuals * 2 samples per individual = 3428
nrow(HeatMAP)==nrow(organs) # TRUE

# Add a column (GeneSet) that categorize if the gene in the "gene" column is SLC25 or not
HeatMAP <- HeatMAP %>%
  mutate(GeneSet = if_else(str_detect(gene,"SLC25"),"SLC25","Biomarker"))

# Add the median of the count value (log2(CPM)) for each gene if it was sampled
# from the same person, same organ, and TCA cycle activity score 
CorDF <- HeatMAP %>%
  group_by(sample_ID, organ, activity, GeneSet) %>%
  summarise(AvgExp = median(qnorm_voom_logCPM))

# Correlation of the SLC25 genes against biomarkers

# Get the names of all the SLC25 genes
# SLCgenes<-unique(HeatMAP$gene[which(HeatMAP$GeneSet=="SLC25")])

# Subset the biomarker genes, group by sample ID, and then take the median count value
# to get the median of the biomarker gene expression 
HeatMAP %>% filter(GeneSet == "Biomarker") %>%
  group_by(sample_ID) %>%
  summarise(medBio = median(qnorm_voom_logCPM)) -> MedianBioMarkerExpressionTissue

# Iterate through the list of SLC25 genes, subset the SLC25 genes, group by the sample ID, 
# combine with the organ value the gene was samples from, then calculate the
# correlation between the median biomarker value between heart and liver
lapply(SLC, function(x){
    HeatMAP %>%
    filter(gene == x) %>%
    left_join(MedianBioMarkerExpressionTissue) %>%
    group_by(organ, gene) %>%
    summarise(cor = cor(medBio, qnorm_voom_logCPM))
}) %>% do.call(rbind,.) -> CorBySLCgene

# Make the dataframe required to make the histogram of correlation values
# of the SLC25 genes
CorBySLCgene %>%
  filter(organ =="heart") %>%
  arrange(cor) %>%
  pull(gene) -> GenePlotOrder

# Barplots of SLC25 correlation values
# Plot the correlation for between each SLC25 gene and the median biomarker
# expression in that organ next to each other; color by organ
CorBySLCgene %>%
  mutate(gene = factor(gene, levels = GenePlotOrder)) %>%
  ggplot(aes(x = gene, y = cor, fill = organ)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = rev(c("#00ACED", "#ED4100"))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "correlation (R) with metabolic biomarkers", x = "SLC25 gene")

# Plot the correlation between SLC25 genes and median metabolic biomarkers
# Plot the correlation values in the heart as the top plot 
# and liver as the bottom
CorBySLCgene %>%
  mutate(gene = factor(gene, levels = GenePlotOrder)) %>%
  ggplot(aes(x = gene, y = cor, fill = organ)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = rev(c("#00ACED", "#ED4100"))) +
  facet_wrap(~organ, nrow = 2) + # this organizes the plots vertically
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  labs(y = "correlation (R) with metabolic biomarkers", x = "SLC25 gene")

# Plot the correlation values (Rs) as side-by-side plots
# PLot on the left is the correlations in heart tissue and on the right is liver
PlotA <- CorBySLCgene %>%
  filter(organ == "heart") %>%
  arrange(cor) %>%
  mutate(gene = factor(gene, levels = unique(gene))) %>%
  ggplot(aes(x = gene, y = cor, fill = organ)) +
  geom_bar(stat = "identity", position = position_dodge(), fill = "#ED4100") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "correlation (R) with metabolic biomarkers", x = "SLC25 gene") +
  ggtitle("heart")

PlotB <- CorBySLCgene %>%
  filter(organ == "liver") %>%
  arrange(cor) %>%
  mutate(gene = factor(gene, levels = unique(gene))) %>%
  ggplot(aes(x = gene, y = cor, fill = organ)) +
  geom_bar(stat = "identity", position = position_dodge(), fill = "#00ACED") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "correlation (R) with metabolic biomarkers", x = "SLC25 gene") +
  ggtitle("liver")

PlotA+PlotB&theme(legend.position = "none") # This plots the barplot side-by-side

# Test
CorBySLCgene%>%
  ggplot(aes(x=organ,y=cor,fill=organ))+
  geom_violin()+
  theme_classic()

# each gene with activity score
AllGenes <- unique(HeatMAP$gene)

# Calculate the correlation between TCA cycle activity and the qnorm + log2CPM
# transformed counts for every gene
# Iterate through the total list of genes (SLC25 + biomarkers) and calculate
# the correlation with TCA cycle activity after grouping by organ
lapply(AllGenes, function(x){
    HeatMAP %>%
    filter(gene == x) %>%
    group_by(organ, GeneSet) %>%
    summarise(cor = cor(activity, qnorm_voom_logCPM)) %>%
    mutate(gene = x)
}) %>% do.call(rbind,.) -> CorActivitybygeneByTissue

# Check the number of expected correlation values are correct
class(CorActivitybygeneByTissue)
nrow(CorActivitybygeneByTissue)==nrow(ac) # 234

# Make a list of the order of the gene names in the heart samples to 
# create the order of the genes
CorActivitybygeneByTissue %>%
  filter(organ == "heart") %>%
  arrange(cor) %>%
  pull(gene) -> GenePlotOrderActivity

# Barplot of the correlation (R) values for the biomarkers and SLC25 next to each other
# color the bars by organ
CorActivitybygeneByTissue %>%
  # mutate() creates new cols that are functions of existing variables
  mutate(gene = factor(gene, levels = GenePlotOrderActivity)) %>% 
  ggplot(aes(x = gene, y = cor, fill = organ)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() +
  ylab("correlation (R) with TCA cycle activity score") +
  ggtitle("Correlation of metabolic biomarkers and SLC25 with TCA cycle activity score") +
  theme(axis.text.x = element_text(angle=90)) +
  facet_wrap(~GeneSet, scales = "free_x")


# Scatter plots with the correlation and a line
# median SLC25 expression in heart vs TCA cycle activity score
M_scatter_heart <- ggpubr::ggscatter(data = CorDF %>% filter(GeneSet=="SLC25" & organ=="heart"), x = "AvgExp", y="activity", add = "reg.line") +
  ggpubr::stat_cor(method = "spearman", label.x = 0.45, label.y=0.9) +
  ggtitle("heart", subtitle = "SLC genes") +
  xlab("median gene set expression") +
  ylab("activity score")

# median biomarker expression vs TCA cycle activity score in heart
O_scatter_heart <- ggpubr::ggscatter(data = CorDF %>% filter(GeneSet=="Biomarker" & organ=="heart"), x = "AvgExp", y="activity", add = "reg.line") +
  ggpubr::stat_cor(method = "spearman", label.x = 0.45, label.y=0.9) +
  ggtitle("",subtitle = "Known biomarker genes") +
  xlab("median gene set expression") +
  ylab("activity score")

# median SLC25 expression in liver vs TCA cycle activity score
M_scatter_Liver <- ggpubr::ggscatter(data = CorDF %>% filter(GeneSet=="SLC25" & organ=="liver"), x = "AvgExp", y="activity", add = "reg.line") +
  ggpubr::stat_cor(method = "spearman", label.x = 0.45, label.y=0.9) +
  ggtitle( "liver",subtitle = "SLC genes") +
  xlab("median gene set expression") +
  ylab("activity score")

# median biomarker expression vs TCA cycle activity score in liver
O_scatter_Liver <- ggpubr::ggscatter(data = CorDF %>% filter(GeneSet=="Biomarker" & organ=="liver"), x = "AvgExp", y="activity", add = "reg.line")+
  ggpubr::stat_cor(method = "spearman", label.x = 0.45, label.y=0.9) +
  ggtitle("",subtitle = "Known biomarker genes") +
  xlab("median gene set expression") +
  ylab("activity score")

# Plot the correlations in the following order
{M_scatter_heart+O_scatter_heart}/{M_scatter_Liver+O_scatter_Liver}

# Create a dataframe with the median expression of SLC25 expression (M_exp) and biomarker expression (O_exp)
CompCor <- CorDF %>% filter(GeneSet == "SLC25") %>% dplyr::select(sample_ID, organ, activity, M_exp = AvgExp) %>%
  left_join(CorDF %>% filter(GeneSet == "Biomarker") %>% dplyr::select(sample_ID, organ, activity, O_exp = AvgExp))

# Correlation between SLC25 and known biomarkers
Exp_scatter<-ggpubr::ggscatter(data = CompCor, x = "M_exp", y ="O_exp", add = "reg.line") +
  ggpubr::stat_cor(method = "spearman") +
  xlab("SLC25 gene set") +
  ylab("Known biomarkers")
Exp_scatter

# Correlation between the average expression of the biomarkers/SLC25 and activity
ggplot(data = CorDF, aes(x = AvgExp, y = activity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(method = "spearman") +
  facet_wrap(~GeneSet)

# Add a column with the organ to the activity dataframe
activity2 <- merge(activity, organs[ ,c("sample_ID", "organ")], by="sample_ID") %>% unique()
activity2 <- unique(activity2)
activity3 <- activity2 %>% filter(str_detect(pathway,"tca"))

# Violin plot of the average SLC25 + biomarker expression per samples vs TCA cycle activity score 
# Violin plot shows a lot of TCA variation in heart relative to liver
activity3 %>%
  dplyr::select(organ, sample_ID, activity) %>%
  unique() %>%
  ggplot(aes(x = organ,y = activity,fill = organ)) +
  geom_violin() +
  geom_jitter(alpha=0.6,width=.1) +
  theme_classic() +
  ylab("TCA cycle activity score") +
  labs(x = "organ, N=148 samples per organ") +
  ggtitle("Activity score by patient and tissue")

# Calculate the mean, standard deviation and coefficient of variation
# increased CoV in heart 
activity3 %>%
  group_by(organ, sample_ID) %>%
  summarise(AvgAct = mean(activity)) %>%
  ungroup() %>%
  group_by(organ) %>%
  summarise(m = mean(AvgAct), sd = sd(AvgAct)) %>%
  mutate(CoV = sd/m)

# Calculate the difference in activity score between heart and liver for each individual
H_act <- HeatMAP %>% filter(organ=="heart") %>% dplyr::select(sample_ID, H_act = activity) %>% unique()  
L_act <- HeatMAP %>% filter(organ == "liver") %>% dplyr::select(sample_ID, L_act = activity) %>% unique()
mutate(diff = H_act[["H_act"]]-L_act[["L_act"]]) %>% arrange(diff)

# Summarize the frequency of the activity score for each organ
table(HeatMAP$organ)

# Make a heatmap of the high, medium, and low expression values
HeatMAP %>%
  #filter(Strata != "med") %>%
  ggplot(aes(x = gene, y = sample_ID, fill = qnorm_voom_logCPM)) +
  geom_tile() +
  facet_grid(organ ~ Strata, scales = "free_y") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text = element_text(size = 3))

# Make a heatmap of the high, medium, and low expression values
HeatMAP %>%
  #filter(Strata != "med") %>%
  group_by(organ=="heart") %>%
  ggplot(aes(x = gene, y = sample_ID, fill = qnorm_voom_logCPM)) +
  geom_tile() +
  facet_grid(gene=="Biomarker" ~ Strata, scales = "free_y") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text = element_text(size = 3))

# Make a heatmap of the high, medium, and low expression values
HeatMAP %>%
  #filter(Strata != "med") %>%
  group_by(organ=="liver") %>%
  ggplot(aes(x = gene, y = sample_ID, fill = qnorm_voom_logCPM)) +
  geom_tile() +
  facet_grid(gene=="SLC25" ~ Strata, scales = "free_y") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text = element_text(size = 3))

# Double check that plotted only the TCA cycle activiy scores
activity$pathway[str_which(activity$pathway,"tca")]
any(is.na(HeatMAP))==TRUE # FALSE

# Make a heatmap of the high, medium, and low expression values
HeatMAP %>%
  #filter(Strata != "med") %>%
  ggplot(aes(x = gene, y = sample_ID, fill = qnorm_voom_logCPM)) +
  geom_tile() +
  facet_grid(organ ~ Strata, scales = "free_y") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text = element_text(size = 3))


# Calculate correlation between biomarker and slc25 gene
# Make a heatmap of the correlation between SLC25 genes and metabolic biomarkers
SLCgene<-unique(HeatMAP$gene[str_which(HeatMAP$gene,"SLC25")])
SLCgene<-SLCgene[!SLCgene=="SLC25A6"]
Biomarkergene<-unique(HeatMAP$gene[which(HeatMAP$GeneSet=="Biomarker")])

# Calculate the correlation for every combination of biomarker and SLC25 gene
# Produces a dataframe
BiomarkerBygeneCor<-lapply(Biomarkergene,function(BG){
  lapply(SLCgene,function(SG){
    cat(BG,SG,'\n')
    HeatMAP%>%
      filter(gene%in%c(BG,SG))%>%
      reshape2::dcast(sample_ID+organ~gene,value.var = "qnorm_voom_logCPM")%>%
      setNames(c("sample_ID","organ","G1","G2"))%>%
      group_by(organ)%>%
      summarise(Cor=cor(G1,G2),G1=BG,G2=SG)
        })%>%do.call(rbind,.)
})%>%do.call(rbind,.)

# Make a heatmap of the biomarker x SLC25 gene correlation values
BiomarkerBygeneCor %>%
  ggplot(aes(G1,G2,fill=Cor)) +
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red")+
  ylab("SLC25 genes") +
  xlab("Biomarker genes") +
  ggtitle("Correlation of SLC25 genes with metabolic biomarkers")

# Calculate the correlation between the TCA activity and each gene
ConditionTest <- reshape2::dcast(data = HeatMAP, formula = activity ~ gene, value.var = "qnorm_voom_logCPM")

# Drop the extra "-" at the end of the gene names
#colnames(ConditionTest) <- str_remove(colnames(ConditionTest), "-")

# Make a list of the covariates (gene names)
covars <- colnames(ConditionTest)[-1]

# Make a dataframe of the effect size, p-value, and SLC25 gene tested vs TCA cycle activity score
sel_colvars <- ""
SLC25_DF <- do.call(rbind,lapply(SLC,function(C){
  cat(C,'\n')
  #SelC <- str_remove(paste(C, sel_colvars, sep = "+"), "\\+$")
  mod <- lm(formula = paste0("activity~", C), data = ConditionTest)
  coef(summary(mod))[,c(1,4)] %>% `colnames<-` (c("estimate", "pvalue")) %>% data.frame() %>% mutate(Test = row.names(.)) %>%
    filter(str_detect(Test,C))
}))

# Change the name of "MT-ND1" to "MTDN1"
colnames(ConditionTest)[41] <- c("MTND1")
colnames(ConditionTest)[41]

# Make a dataframe of the effect size, p-value, and biomarker gene tested vs TCA cycle activity score
Biomarker_DF <- do.call(rbind,lapply(biomarkers,function(C){
  cat(C,'\n')
 # SelC <- str_remove(paste(C, sel_colvars, sep = "+"), "\\+$")
  mod <- lm(formula = paste0("activity~", C), data = ConditionTest)
  coef(summary(mod))[,c(1,4)] %>% `colnames<-` (c("estimate", "pvalue")) %>% data.frame() %>% mutate(Test = row.names(.)) %>%
    filter(str_detect(Test,C))
}))

# Combine the tables together
combined_DF <- rbind(Biomarker_DF, SLC25_DF)
row.names(combined_DF) <- NULL
combined_DF <- combined_DF %>% relocate("Test")
head(combined_DF)
tail(combined_DF)

# Write to file
write.table(combined_DF, "~/Downloads/correlation_SLC25_biomarkers_TCA_activity_score.tsv", sep="\t", row.names = FALSE)

# Histogram of the frequency of TCA cycle activity in heart and liver
HeatMAP%>%
  dplyr::select(sample_ID,organ,activity)%>%
  unique()%>%
  ggplot(aes(x=activity,fill=organ))+
  geom_histogram(binwidth = .025,color="black")+
  facet_wrap(~organ,ncol=1)+
  theme_classic()+
  labs(x="TCA cycle activity score \n N=148 individuals")+
  ggtitle("Histogram of TCA cycle activity scores in heart and liver")

