#R file for the transcriptomic analysis of Dae D. Chung's paper: 
#"Sex Differences in the Transcriptome of Extracellular Vesicles Secreted by Fetal Neural Stem Cells and Effects of Chronic Alcohol Exposure"
#Details for the RNA Isolation, RNA-seq Library Preparation, and Sequencing are written in the Materials and Methods section.
#Raw (fastq) and processed data (htseq) files for all samples' total RNA are deposited in NCBI repository under accession number GSE214545.


#Before anything, get where your directory is and set your working directory. 
getwd()
setwd()

#side note: Deseq2 analysis comparing EV vs. cell samples and female vs. male samples were both done on
#the Galaxy platform, and not in R. So I do not have R codes for those. 
#We will be importing Deseq2 analyzed files from Galaxy to R.

#First, we will import deseq2 data that did differential expression analysis, comparing EV samples to cell samples.
#EV vs Cell first####
#import dataset from galaxy
#these are output dataset from deseq2 done from Galaxy platform, where the input data is from htseq. 
#the two files i get are the results file and the norm_counts file
#the result file has different columns with
#Column	Description
#1	    Gene Identifiers
#2	    mean normalised counts, averaged over all samples from both conditions
#3	    the logarithm (to basis 2) of the fold change (See the note in inputs section)
#4	    standard error estimate for the log2 fold change estimate
#5	    Wald statistic
#6	    p value for the statistical significance of this change
#7	    p value adjusted for multiple testing with the Benjamini-Hochberg procedure which controls false discovery rate (FDR)

#mean normalized counts are same for files that use all 36 samples, but will be different if only some of the 36 samples are used,
#like 120 vs 0, or female ev vs male ev, etc.

#norm_count files with normalized counts, given by DESeq2 (built-in normalization method)
#https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#The normalization method is by DESeq2's median of ratios, which is 
#counts divided by sample-specific size factors 
#determined by median ratio of gene counts relative to geometric mean per gene.
#This accounts for sequencing depth (necessary for comparison of gene expression between samples, like if sample A has doubled in expression for most genes vs Sample B) 
#and RNA composition (recommended for accurate comparison of expression between samples, and is particularly important when performing differential expression analyses).
#This normalization is recommended for gene count comparisons between samples and for DE analysis; NOT for within sample comparisons,
#because it does not account for different gene length of genes within samples.

#before importing to R, I decided to put in pregnancy, EV_vs_Cell, and Treatment rows in excel.

#row names are currently numeric, need to make rownames = Gene_ID
#need to make it a data frame to work with it
EV_vs_Cell_results<-as.data.frame(DESeq2_result_file_EV_vs_Cell_RNA_Star_gene_Galaxy1642)
EV_vs_Cell_norm_counts<-as.data.frame(DESeq2_normalized_counts_file_EV_vs_Cell_RNA_Star_gene_Galaxy1644_Normalized_counts_file_on_data_1638_data_1632_and_others_)
View(EV_vs_Cell_results)
View(EV_vs_Cell_norm_counts)
#make first column into row names
row.names(EV_vs_Cell_results) <- EV_vs_Cell_results[,1]
row.names(EV_vs_Cell_norm_counts) <- EV_vs_Cell_norm_counts[,1]
#delete first column
EV_vs_Cell_results[1] <- NULL
EV_vs_Cell_norm_counts[1] <- NULL



#transpose so that columns are now Gene_ID and rows are samples
#and need to make it a data frame to work with it
tEV_vs_Cell_norm_counts<-as.data.frame(t(EV_vs_Cell_norm_counts))
View(tEV_vs_Cell_norm_counts)
#convert all columns of normalized counts to numeric from factor, except first 4 columns that describe sample traits
tEV_vs_Cell_norm_counts[,5:ncol(tEV_vs_Cell_norm_counts)] <- 
  lapply(tEV_vs_Cell_norm_counts[,5:ncol(tEV_vs_Cell_norm_counts)], function(x) as.numeric(as.character(x)))
#get rid of columns that have no variance ans so can be used in prcomp
which(apply(tEV_vs_Cell_norm_counts, 2, var)==0)
clean_tEV_vs_Cell_norm_counts<-tEV_vs_Cell_norm_counts[ , which(apply(tEV_vs_Cell_norm_counts, 2, var) != 0)]
#convert all columns of normalized counts to numeric from factor
clean_tEV_vs_Cell_norm_counts[,] <- lapply(clean_tEV_vs_Cell_norm_counts[,], function(x) as.numeric(as.character(x)))
View(clean_tEV_vs_Cell_norm_counts[,1:30])
#get ride of columns with sample traits
clean_tEV_vs_Cell_norm_counts<-clean_tEV_vs_Cell_norm_counts[,-c(1:2)]
View(clean_tEV_vs_Cell_norm_counts[,1:30])

####PCA######
install.packages("ggplot2")
install.packages("grid")
install.packages("gridExtra")
install.packages("ggfortify")
library(ggplot2)
library(grid)
library(gridExtra)
library(ggfortify)

#PCA of top 500 most variant genes#
#https://www.datacamp.com/tutorial/pca-analysis-r

install.packages("matrixStats")
library(matrixStats)
#need it to have no NAs and genes are rows and samples are columns
var_EV_vs_Cell <- as.data.frame(t(clean_tEV_vs_Cell_norm_counts))
View(var_EV_vs_Cell[,1:36])

#take row var (var across gene) and add as column at end
var_EV_vs_Cell$row_var <- rowVars(as.matrix(var_EV_vs_Cell[,c(1:36)]))
#sort table by variance, negative sign so largest to smallest, then select highest 500
sorted_var_EV_vs_Cell <- var_EV_vs_Cell[order(-var_EV_vs_Cell$row_var),][1:500,]
View(sorted_var_EV_vs_Cell)
#transpose back to rows subjects and columns genes for PCA and get rid of row_var
t_var_EV_vs_Cell <- as.data.frame(t(sorted_var_EV_vs_Cell[,-37]))
PCA_var_EV_vs_Cell<-prcomp(t_var_EV_vs_Cell, center = TRUE,scale. = TRUE)
summary(PCA_var_EV_vs_Cell)
View(PCA_var_EV_vs_Cell)


ggplot(PCA_var_EV_vs_Cell)

dtp <- data.frame('EV_vs_Cell' = tEV_vs_Cell_norm_counts$EV_vs_Cell, 
                  'Sex' = tEV_vs_Cell_norm_counts$Sex,
                  'Pregnancy' = tEV_vs_Cell_norm_counts$Pregnancy,
                  'Treatment' = tEV_vs_Cell_norm_counts$Treatment,
                  PCA_var_EV_vs_Cell$x[,1:2]) # the first two components are selected (NB: you can also select 3 for 3D plottings or 3+)
View(dtp)

#by EV_vs_Cell
png("EV_vs_Cell_PCA_top_500_variant_by_EV_vs_Cell.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp) + 
  geom_point(aes(x = PC1, y = PC2, col = Sex), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = EV_vs_Cell)) +
  scale_color_manual(values = c("#D55E00", "#000000", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00"))
dev.off()
getwd()

autoplot(PCA_var_EV_vs_Cell)
autoplot(PCA_var_EV_vs_Cell, loadings=TRUE, loadings.colour='darkorchid4', loadings.label=TRUE, loadings.label.size=3)


#autoplot shows the percentage of pc1 and pc2, which were 52.98 and 9.26 in this case.
#so graph with those in x and y labels 


png("EV_vs_Cell_PCA_top_500_variant_by_EV_vs_Cell_variable_labels.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp) + 
  geom_point(aes(x = PC1, y = PC2, col = Sex), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = EV_vs_Cell)) +
  xlab("PC 1 (52.98%)") + 
  ylab("PC 2 (9.26%)") +
  scale_color_manual(values = c("#D55E00", "#000000", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00"))+
  scale_y_reverse()+ scale_x_reverse()
dev.off()



#by Sex
png("EV_vs_Cell_PCA_top_500_variant_by_Sex.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp) + 
  geom_point(aes(x = PC1, y = PC2, col = EV_vs_Cell), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = Sex)) +
  scale_color_manual(values = c("#D55E00", "#000000", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00"))
dev.off()

#samples grouped by EV and Cell


PCA_var_EV_vs_Cell_components <- as.data.frame(PCA_var_EV_vs_Cell$rotation)
View(PCA_var_EV_vs_Cell_components)
write.csv(PCA_var_EV_vs_Cell_components, "EV_vs_Cell PCA contibutors.csv")



####Volcano Plot of EV_vs_Cell_results####
#if not done already then:
#https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

install.packages("ggplot2")
install.packages("ggrepel")

library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

EV_vs_Cell_results<-as.data.frame(DESeq2_result_file_EV_vs_Cell_RNA_Star_gene_Galaxy1642)
View(EV_vs_Cell_results)
#make first column into row names
row.names(EV_vs_Cell_results) <- EV_vs_Cell_results[,1]
#delete first column
EV_vs_Cell_results[1] <- NULL
colnames(EV_vs_Cell_results)
View(EV_vs_Cell_results)
View(EV_vs_Cell_results_sig_genes)

write.csv(EV_vs_Cell_results, "DESEQ2_results_EV_vs_Cell.csv")
write.csv(EV_vs_Cell_results_sig_genes, "DESEQ2_results_EV_vs_Cell_sig_genes.csv")


#adj-P significantly differentially expressed genes
length(which(EV_vs_Cell_results$`P-adj`<0.05))
#21539

#cool, lets see how many of them are downregulated (meaning less in EV vs cell, log2(FC)<0) vs upregulated
length(which(EV_vs_Cell_results$`P-adj`<0.05 & EV_vs_Cell_results$`log2(FC)`<0))
#10422
length(which(EV_vs_Cell_results$`P-adj`<0.05 & EV_vs_Cell_results$`log2(FC)`>0))
#11117


#The default cut-off for log2FC is >|1|, if FC is 0.5 then you get -1, if 2 then +1, if 1 then 0; 
#the default cut-off for P value is 10e-6. 
pdf("EV_vs_Cell_Volcano_plot_default.pdf")
EnhancedVolcano(EV_vs_Cell_results,
                lab = rownames(EV_vs_Cell_results),
                title = 'EV vs. Cell',
                x = 'log2(FC)',
                y = 'P-adj',
                xlim = c(-6, 6))
dev.off()

#custom cut-offs: FC_log2 0.5-1.5; default adj_pval of 10e-6
pdf("ANOVA_EV_vs_Cell_Volcano_plot_custom.pdf")
EnhancedVolcano(EV_vs_Cell_results,
                lab = rownames(EV_vs_Cell_results),
                x = 'log2(FC)',
                y = 'P-adj',
                xlim = c(-2, 2),
                title = 'EV vs. Cell',
                pCutoff = 10e-6,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30')
dev.off()

#####PCA components####### 
#so far:
#1. pca shows samples grouped by ev vs cell, with ev samples grouping together on the negative pc1 side/value,
#and cell smaples grouping together on the positive pc1 side/value.
#2. volcano plot of deseq2 results for ev vs cell shows bunch of transcripts that are significantly different,
#3. intersect sig adj pvalue genes from deseq2 and pc1 contributors from the top 500 variants.
#4. focus on the largest contributors of PC1 from these genes, with 
#4a. order by largest negative PC1 values to look at genes that contributed most in EV samples being grouped away from cell samples
#that contributed most to the negative value of avg PC1 for EV samples
#4b. order by largest positive PC1 values to look at genes that contributed most in cell samples being grouped away from EV samples
#that contributed most to the positive value of avg PC1 for Cell samples
#4c pathway analysis



#3. intersect sig adj pvalue genes from deseq2 and pc1 contributors from the top 500 variants.####

View(EV_vs_Cell_results)
#only include transcripts that have adjusted pvalue of 0.05 and less
EV_vs_Cell_results_sig_genes<-as.data.frame(EV_vs_Cell_results[which(EV_vs_Cell_results$`P-adj`<0.05),])
View(EV_vs_Cell_results_sig_genes)

#then order with log2(FC) that are positive, since that means greater in EV samples than cell samples.
#the default of command "order()" is smallest to largest, so negative value to positive value descending.
#thus, we do order(-), to reverse that and get the highest positive values on the top. 
EV_vs_Cell_results_sig_genes <- EV_vs_Cell_results_sig_genes[order(-EV_vs_Cell_results_sig_genes$`log2(FC)`),]
#this d.f. does not have a column with Gene_ID, so make rownames into a column
EV_vs_Cell_results_sig_genes$Gene_ID<-rownames(EV_vs_Cell_results_sig_genes)


PC1_genes <- as.data.frame(PCA_var_EV_vs_Cell_components[order(-PCA_var_EV_vs_Cell_components$PC1),][,])
PC1_genes$Gene_ID <- rownames(PC1_genes)
PC1_sig_genes <- as.data.frame(intersect(PC1_genes$Gene_ID, EV_vs_Cell_results_sig_genes$Gene_ID))
colnames(PC1_sig_genes)[1] <- "Gene_ID"
#add PC1 contribution and Fold change and adj p-value
row.names(PC1_sig_genes) <- PC1_sig_genes$Gene_ID
PC1_sig_genes <- merge(PC1_sig_genes,PCA_var_EV_vs_Cell_components["PC1"],by="row.names",all.x=TRUE)
PC1_sig_genes[1] <- NULL
row.names(PC1_sig_genes) <- PC1_sig_genes$Gene_ID
EV_vs_Cell_results_sig_genes <- as.data.frame(EV_vs_Cell_results_sig_genes)
row.names(EV_vs_Cell_results_sig_genes) <- EV_vs_Cell_results_sig_genes$Gene_ID
PC1_sig_genes <- merge(PC1_sig_genes,EV_vs_Cell_results_sig_genes["log2(FC)"],by="row.names",all.x=TRUE)
PC1_sig_genes[1] <- NULL
row.names(PC1_sig_genes) <- PC1_sig_genes$Gene_ID
PC1_sig_genes <- merge(PC1_sig_genes,EV_vs_Cell_results_sig_genes["P-adj"],by="row.names",all.x=TRUE)
PC1_sig_genes[1] <- NULL

getwd()
setwd("C:/Users/daehy/Dropbox/Transcriptomic/R")
View(PC1_sig_genes)

write.csv(PC1_sig_genes, "EV_vs_Cell_PC1_contributor_sig_genes.csv")

###4c pathway analysis####


####let's do pathway analysis####

#will need to convert to Enterz ID
#https://yulab-smu.github.io/clusterProfiler-book/chapter14.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
#install db for mouse IDs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DOSE")
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(DOSE)

library(ReactomePA)
#create vector of IDs
View(PC1_sig_genes)
#keytypes to list all supporting types
keytypes(org.Mm.eg.db)
ID_location <- PC1_sig_genes$Gene_ID
eg_location <- bitr(ID_location, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg_location)
#  SYMBOL         ENTREZID
#1 1500009L16Rik    69784
#2 2410006H16Rik    69221
#3 4930527F14Rik    67651
#4 5031425E22Rik   269630
#5 A730017L22Rik   613258
#6         Abca1    11303


#merge PC1_sig_genes$log2(FC) and eg_location
#btw, all negative PC1 values have positive log2(FC) and all positive PC1 values have negative log2(FC)



#rename column in eg so matching title
View(eg_location)
colnames(eg_location)[1] <- "Gene_ID"



####let's do pathway analysis for EV_vs_Cell_results_sig_genes, adjusted pvalue of 0.05####
View(EV_vs_Cell_results_sig_genes)
ID_location <- EV_vs_Cell_results_sig_genes$Gene_ID
eg_location <- bitr(ID_location, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg_location)
#  SYMBOL         ENTREZID
#1 1500009L16Rik    69784
#2 2410006H16Rik    69221
#3 4930527F14Rik    67651
#4 5031425E22Rik   269630
#5 A730017L22Rik   613258
#6         Abca1    11303


#merge PC1_sig_genes$log2(FC) and eg_location
#btw, all negative PC1 values have positive log2(FC) and all positive PC1 values have negative log2(FC)



#rename column in eg so matching title
View(eg_location)
colnames(eg_location)[1] <- "Gene_ID"

#merge based on Gene_ID for log2(FC)
EV_vs_Cell_Entrez_log2_FC <- 
  merge(eg_location, EV_vs_Cell_results_sig_genes[,c(2,7)], by = "Gene_ID")
View(EV_vs_Cell_Entrez_log2_FC)
#get rid of gene id so i have entrez and fc for 2 columns only
EV_vs_Cell_Entrez_log2_FC <- EV_vs_Cell_Entrez_log2_FC[,-1]
#make it numeric like geneList
#https://bioinformatics.stackexchange.com/questions/13036/reactomepa-dataset
myGeneList <- EV_vs_Cell_Entrez_log2_FC[,2]
names(myGeneList) <- as.character(EV_vs_Cell_Entrez_log2_FC[,1])
myGeneList <- sort(myGeneList, decreasing = TRUE)
head(myGeneList)
View(myGeneList)

#start processing for Reactome
#limiting genes with log2 fold change that are greater than 0.
#If we use log2(fold change), fold changes lower than 1 (when B > A) become negative, 
#while those greater than 1 (A > B) become positive. 
#Now the values are symmetrical and it's easier to see fold changes in both directions on one plot.
#so, we will do anything above 0 as enriched in EVs since that is fold change of greater than 1. 
#use abs() if you dont care if the fc is negative or positive.

#data(geneList, package="DOSE")
#View(geneList)
class(myGeneList) #its numeric
de <- names(myGeneList)[myGeneList > 0]
head(de)

#The enrichPathway function allows users to select an appropriate background of genes as the baseline. 
#The gsePathway function supports GSEA to evaluate the enriched reactome pathways of high-throughput data. 
#These approaches can be used to verify interesting altered pathways and to identify unanticipated pathway associations.
x <- enrichPathway(gene=de, organism = "mouse", pvalueCutoff=0.05, pAdjustMethod = "BH", readable=T)
head(as.data.frame(x))
er <- as.data.frame(x)
#enriched term found. 

getwd()
write.csv(er, "DESEQ2_results_EVenriched_EV_vs_Cell_foldchange2_Reactome_Results.csv")

png("DESEQ2_results_EVenriched_EV_vs_Cell_foldchange2_bar_Reactome.png", width = 4000, height = 2000, res = 300)
barplot(x, showCategory=10)
dev.off()

png("DESEQ2_results_EVenriched_EV_vs_Cell_foldchange2_dotplot_Reactome.png", width = 2000, height = 2000, res = 300)
dotplot(x, showCategory=10)
dev.off()

png("DESEQ2_results_EVenriched_EV_vs_Cell_foldchange2_dotplot_Reactome.png", width = 1500, height = 2500, res = 300)
dotplot(x, showCategory=15, font.size = 10)
dev.off()



##https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
##additional dot plot adjustments!!!!!!

x2<- pairwise_termsim(x) 
emapplot(x2, color="pvalue")

png("DESEQ2_results_EVenriched_EV_vs_Cell_foldchange2_enrichMap_Reactome.png", width = 2000, height = 2000, res = 300)
emapplot(x2)
dev.off()

png("DESEQ2_results_EVenriched_EV_vs_Cell_foldchange2_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         cex_label_gene = 0.8)
dev.off()


png("DESEQ2_results_EVenriched_EV_vs_Cell_foldchange2_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         node_label="category")
dev.off()

#GSEA
#9.3 Reactome pathway gene set enrichment analysis; gsePathway()
#Gene Set Enrichment Analysis (GSEA) is a computational method 
#that determines whether an a priori defined set of genes shows statistically
#significant, concordant differences between two biological states (e.g. phenotypes).
y <- gsePathway(myGeneList, 
                organism = "mouse",
                pvalueCutoff=0.05,
                pAdjustMethod="BH", verbose=FALSE)
y2 <- pairwise_termsim(y) 

png("DESEQ2_results_EV_vs_Cell_foldchange2_GSEA_Reactome.png", width = 3000, height = 2000, res = 300)
emapplot(y2)
dev.off()

#gsePathway works, because it has genes enriched in EVs or cells in the "myGeneList". 
#whenever enrichment score and NES scores are negative, 
#that means they are genes that are at the bottom of the list or those with negative log2(FC) scores for my case.
#this means i have results for over-represented pathways for EV enriched genes and for cell enriched genes. 

head(as.data.frame(y2))
View(as.data.frame(y))
head(as.data.frame(y2))
er <- as.data.frame(y2)

getwd()
write.csv(er, "DESEQ2_results_EVenriched_EV_vs_Cell_foldchange2_GSEA_Reactome.csv")

#bar and dotplots are useless for looking at ev enriched genes unless i set the order by enrichment score,
#where positive is ev enriched.
png("DESEQ2_results_EVenriched_EV_vs_Cell_foldchange2_bar_GSEA_Reactome.png", width = 4000, height = 2000, res = 300)
barplot(y2, showCategory=10)
dev.off()


png("DESEQ2_results_EVenriched_EV_vs_Cell_foldchange2_dotplot_GSEA_Reactome.png", width = 1500, height = 2500, res = 300)
dotplot(y2, showCategory=15, font.size = 10)
dev.off()



##let's try cell enriched:####
de_cell <- names(myGeneList)[myGeneList < 0]

head(de_cell)
View(de_cell)
#de_cell<-rev(de_cell) do not need to reverse row order; does not affect the pathways analysis

#The enrichPathway function allows users to select an appropriate background of genes as the baseline. 
#The gsePathway function supports GSEA to evaluate the enriched reactome pathways of high-throughput data. 
#These approaches can be used to verify interesting altered pathways and to identify unanticipated pathway associations.
x <- enrichPathway(gene=de_cell, organism = "mouse", pvalueCutoff=0.05, pAdjustMethod = "BH", readable=T)
head(as.data.frame(x))
er <- as.data.frame(x)
#enriched term found for cells. 

getwd()
write.csv(er, "DESEQ2_results_Cellenriched_EV_vs_Cell_foldchange2_Reactome_Results.csv")

png("DESEQ2_results_Cellenriched_EV_vs_Cell_foldchange2_barplot_Reactome.png", width = 4000, height = 2000, res = 300)
barplot(x, showCategory=10)
dev.off()

png("DESEQ2_results_Cellenriched_EV_vs_Cell_foldchange2_dotplot_Reactome.png", width = 2000, height = 2000, res = 300)
dotplot(x, showCategory=10)
dev.off()

png("DESEQ2_results_Cellenriched_EV_vs_Cell_foldchange2_dotplot_Reactome.png", width = 1500, height = 2500, res = 300)
dotplot(x, showCategory=15, font.size = 10)
dev.off()



##https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
##additional dot plot adjustments!!!!!!

x2 <- pairwise_termsim(x) 
emapplot(x2, color="pvalue")

png("DESEQ2_results_Cellenriched_EV_vs_Cell_foldchange2_enrichMap_Reactome.png", width = 2000, height = 2000, res = 300)
emapplot(x2)
dev.off()

png("DESEQ2_results_Cellenriched_EV_vs_Cell_foldchange2_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         cex_label_gene = 0.8)
dev.off()

#png("EV_vs_Cell_Cellenriched_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
#cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
#         cex_label_category = 0.8, 
#         node_label="category")
#dev.off()

#GSEA


#KEGG enrichment analysis
#visualize KEGG

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("pathview")
#library("pathview")

#7.2 KEGG pathway over-representation analysis


View(de_cell)
kk <- enrichKEGG(gene         = de_cell,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
head(kk)

View(as.data.frame(kk))
er_kk<- as.data.frame(kk)
getwd()
write.csv(er_kk, "DESEQ2_results_Cellenriched_EV_vs_Cell_foldchange2_KEGG_Results.csv")



###############sex, female vs male, separated by EV and cell samples#####

#import dataset from galaxy
#these are output dataset from deseq2, where the input data is from htseq. 
#the two files i get are the results file and the norm_counts file
#the result file has different columns with
#Column	Description
#1	    Gene Identifiers
#2	    mean normalised counts, averaged over all samples from both conditions
#3	    the logarithm (to basis 2) of the fold change (See the note in inputs section)
#4	    standard error estimate for the log2 fold change estimate
#5	    Wald statistic
#6	    p value for the statistical significance of this change
#7	    p value adjusted for multiple testing with the Benjamini-Hochberg procedure which controls false discovery rate (FDR)

#mean normalized counts are same for files that use all 36 samples, but will be different if only some of the 36 samples are used,
#like 120 vs 0, or female ev vs male ev, etc.

#norm_count files with normalized counts, given by DESeq2 (built-in normalization method)
#https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#The normalization method is by DESeq2's median of ratios, which is 
#counts divided by sample-specific size factors 
#determined by median ratio of gene counts relative to geometric mean per gene.
#This accounts for sequencing depth (necessary for comparison of gene expression between samples, like if sample A has doubled in expression for most genes vs Sample B) 
#and RNA composition (recommended for accurate comparison of expression between samples, and is particularly important when performing differential expression analyses).
#This normalization is recommended for gene count comparisons between samples and for DE analysis; NOT for within sample comparisons,
#because it does not account for different gene length of genes within samples.


#before importing to R, put in pregnancy, EV_vs_Cell, and Treatment rows in excel.
View(DESeq2_normalized_counts_file_Sex_EV_Cell_RNA_Star_gene_Galaxy1650)
View(DESeq2_result_file_Sex_EV_Cell_RNA_Star_gene_Galaxy1648)
View(DESeq2_normalized_counts_file_Sex_EV_only_RNA_Star_gene_Galaxy1662)
View(DESeq2_result_file_Sex_EV_only_RNA_Star_gene_Galaxy1660)
View(DESeq2_normalized_counts_Sex_Cell_only_RNA_Star_gene_Galaxy1668)
View(DESeq2_result_file_Sex_Cell_only_RNA_Star_gene_Galaxy1666)


#row names are currently numeric, need to make rownames = Gene_ID
#need to make it a data frame to work with it
Sex_EV_Cell_results<-as.data.frame(DESeq2_result_file_Sex_EV_Cell_RNA_Star_gene_Galaxy1648)
Sex_EV_Cell_norm_counts<-as.data.frame(DESeq2_normalized_counts_file_Sex_EV_Cell_RNA_Star_gene_Galaxy1650)
Sex_EV_results<-as.data.frame(DESeq2_result_file_Sex_EV_only_RNA_Star_gene_Galaxy1660)
Sex_EV_norm_counts<-as.data.frame(DESeq2_normalized_counts_file_Sex_EV_only_RNA_Star_gene_Galaxy1662)
Sex_Cell_results<-as.data.frame(DESeq2_result_file_Sex_Cell_only_RNA_Star_gene_Galaxy1666)
Sex_Cell_norm_counts<-as.data.frame(DESeq2_normalized_counts_Sex_Cell_only_RNA_Star_gene_Galaxy1668)
View(Sex_EV_Cell_results)
View(Sex_EV_Cell_norm_counts)
View(Sex_EV_results)
#make first column into row names
row.names(Sex_EV_Cell_results) <- Sex_EV_Cell_results[,1]
row.names(Sex_EV_Cell_norm_counts) <- Sex_EV_Cell_norm_counts[,1]
row.names(Sex_EV_results) <- Sex_EV_results[,1]
row.names(Sex_EV_norm_counts) <- Sex_EV_norm_counts[,1]
row.names(Sex_Cell_results) <- Sex_Cell_results[,1]
row.names(Sex_Cell_norm_counts) <- Sex_Cell_norm_counts[,1]

#delete first column
Sex_EV_Cell_results[1] <- NULL
Sex_EV_Cell_norm_counts[1] <- NULL
Sex_EV_results[1] <- NULL
Sex_EV_norm_counts[1] <- NULL
Sex_Cell_results[1] <- NULL
Sex_Cell_norm_counts[1] <- NULL

View(Sex_Cell_norm_counts)

View(EV_vs_Cell_results)
View(EV_vs_Cell_results_sig_genes)

View(Sex_EV_results)
View(Sex_Cell_results)

write.csv(EV_vs_Cell_results, "DESEQ2_results_EV_vs_Cell.csv")
write.csv(EV_vs_Cell_results_sig_genes, "DESEQ2_results_EV_vs_Cell_sig_genes.csv")

View(Sex_EV_Cell_results)
Sex_EV_Cell_results_sig_genes<-Sex_EV_Cell_results[which(Sex_EV_Cell_results$`P-adj`<0.05),]
View(Sex_EV_Cell_results_sig_genes)
write.csv(Sex_EV_Cell_results_sig_genes, "DESEQ2_results_Sex_EV_and_Cell_sig_genes.csv")

Sex_EV_results_sig_genes<-Sex_EV_results[which(Sex_EV_results$`P-adj`<0.05),]
View(Sex_EV_results_sig_genes)
write.csv(Sex_EV_results, "DESEQ2_results_EV.csv")
write.csv(Sex_EV_results_sig_genes, "DESEQ2_results_EV_sig_genes.csv")

Sex_Cell_results_sig_genes<-Sex_Cell_results[which(Sex_Cell_results$`P-adj`<0.05),]
View(Sex_Cell_results_sig_genes)
write.csv(Sex_Cell_results, "DESEQ2_results_Cell.csv")
write.csv(Sex_Cell_results_sig_genes, "DESEQ2_results_Cell_sig_genes.csv")

#for EV
#adj-P significantly differentially expressed genes
length(which(Sex_EV_results$`P-adj`<0.05))
#2859
#cool, lets see how many of them are downregulated (meaning less in female vs male, log2(FC)<0) vs upregulated
length(which(Sex_EV_results$`P-adj`<0.05 & Sex_EV_results$`log2(FC)`<0))
#1937
length(which(Sex_EV_results$`P-adj`<0.05 & Sex_EV_results$`log2(FC)`>0))
#922

#for cell
length(which(Sex_Cell_results$`P-adj`<0.05))
#2437
#cool, lets see how many of them are downregulated (meaning less in female vs male, log2(FC)<0) vs upregulated
length(which(Sex_Cell_results$`P-adj`<0.05 & Sex_Cell_results$`log2(FC)`<0))
#1638
length(which(Sex_Cell_results$`P-adj`<0.05 & Sex_Cell_results$`log2(FC)`>0))
#799

Sex_EV_results_sig_genes
Sex_Cell_results_sig_genes
View(Sex_EV_results_sig_genes)
View(Sex_Cell_results_sig_genes)
length(which(as.numeric(Sex_EV_results_sig_genes$Gene_ID) & as.numeric(Sex_Cell_results_sig_genes$Gene_ID)))
#merge 2df to see which ones are present in both based on Gene_ID
Sex_EV_results_sig_genes_in_Cell_too <- 
  merge(Sex_EV_results_sig_genes, Sex_Cell_results_sig_genes, by = "Gene_ID", no.dups = TRUE)
View(Sex_EV_results_sig_genes_in_Cell_too)

#you get duplicate columns, so solve that by choosing Genes in Sex_EV_results_sig_genes that are also present in Sex_Cell_results_sig_genes:
Sex_EV_results_sig_genes_in_Cell_too<-subset(Sex_EV_results_sig_genes,Sex_EV_results_sig_genes$Gene_ID %in%Sex_Cell_results_sig_genes$Gene_ID)
View(Sex_EV_results_sig_genes_in_Cell_too)
write.csv(Sex_EV_results_sig_genes_in_Cell_too, "DESEQ2_Sex_EV_results_sig_genes_in_Cell_too.csv")


getwd()
setwd("C:/Users/daehy/Dropbox/Transcriptomic/R")

####PCA######
install.packages("ggplot2")
install.packages("grid")
install.packages("gridExtra")
install.packages("ggfortify")
library(ggplot2)
library(grid)
library(gridExtra)
library(ggfortify)

remove.packages("rlang")
install.packages("rlang")

#transpose so that columns are now Gene_ID and rows are samples
#and need to make it a data frame to work with it
tSex_EV_Cell_norm_counts<-as.data.frame(t(Sex_EV_Cell_norm_counts))
tSex_EV_norm_counts<-as.data.frame(t(Sex_EV_norm_counts))
tSex_Cell_norm_counts<-as.data.frame(t(Sex_Cell_norm_counts))


View(tSex_EV_Cell_norm_counts)
View(tSex_EV_norm_counts)

#convert all columns of normalized counts to numeric from factor, except first 4 columns that describe sample traits
tSex_EV_Cell_norm_counts[,5:ncol(tSex_EV_Cell_norm_counts)] <- 
  lapply(tSex_EV_Cell_norm_counts[,5:ncol(tSex_EV_Cell_norm_counts)], function(x) as.numeric(as.character(x)))
tSex_EV_norm_counts[,5:ncol(tSex_EV_norm_counts)] <- 
  lapply(tSex_EV_norm_counts[,5:ncol(tSex_EV_norm_counts)], function(x) as.numeric(as.character(x)))
tSex_Cell_norm_counts[,5:ncol(tSex_Cell_norm_counts)] <- 
  lapply(tSex_Cell_norm_counts[,5:ncol(tSex_Cell_norm_counts)], function(x) as.numeric(as.character(x)))


#get rid of columns that have no variance ans so can be used in prcomp
which(apply(tSex_EV_Cell_norm_counts, 2, var)==0)
clean_tSex_EV_Cell_norm_counts<-tSex_EV_Cell_norm_counts[ , which(apply(tSex_EV_Cell_norm_counts, 2, var) != 0)]
which(apply(tSex_EV_norm_counts, 2, var)==0)
clean_tSex_EV_norm_counts<-tSex_EV_norm_counts[ , which(apply(tSex_EV_norm_counts, 2, var) != 0)]
which(apply(tSex_Cell_norm_counts, 2, var)==0)
clean_tSex_Cell_norm_counts<-tSex_Cell_norm_counts[ , which(apply(tSex_Cell_norm_counts, 2, var) != 0)]


#convert all columns of normalized counts to numeric from factor
clean_tSex_EV_Cell_norm_counts[,] <- lapply(clean_tSex_EV_Cell_norm_counts[,], function(x) as.numeric(as.character(x)))
clean_tSex_EV_norm_counts[,] <- lapply(clean_tSex_EV_norm_counts[,], function(x) as.numeric(as.character(x)))
clean_tSex_Cell_norm_counts[,] <- lapply(clean_tSex_Cell_norm_counts[,], function(x) as.numeric(as.character(x)))
View(clean_tSex_Cell_norm_counts[,1:30])


#get ride of columns with sample traits
clean_tSex_EV_Cell_norm_counts<-clean_tSex_EV_Cell_norm_counts[,-c(1:2)]
clean_tSex_EV_norm_counts<-clean_tSex_EV_norm_counts[,-c(1:2)]
clean_tSex_Cell_norm_counts<-clean_tSex_Cell_norm_counts[,-c(1:2)]
View(clean_tSex_Cell_norm_counts[,1:30])

#PCA of most variant genes########
install.packages("matrixStats")
library(matrixStats)
#need it to have no NAs and genes are rows and samples are columns
var_Sex_EV_Cell <- as.data.frame(t(clean_tSex_EV_Cell_norm_counts))
var_Sex_EV <- as.data.frame(t(clean_tSex_EV_norm_counts))
var_Sex_Cell <- as.data.frame(t(clean_tSex_Cell_norm_counts))
View(var_Sex_EV_Cell[,1:36])

#take row var (var across gene) and add as column at end
var_Sex_EV_Cell$row_var <- rowVars(as.matrix(var_Sex_EV_Cell[,c(1:36)]))
var_Sex_EV$row_var <- rowVars(as.matrix(var_Sex_EV[,c(1:18)]))
var_Sex_Cell$row_var <- rowVars(as.matrix(var_Sex_Cell[,c(1:18)]))

#sort table by variance, negative sign so largest to smallest, then select highest 500
sorted_var_Sex_EV_Cell <- var_Sex_EV_Cell[order(-var_Sex_EV_Cell$row_var),][1:500,]
sorted_var_Sex_EV <- var_Sex_EV[order(-var_Sex_EV$row_var),][1:500,]
sorted_var_Sex_Cell <- var_Sex_Cell[order(-var_Sex_Cell$row_var),][1:500,]
View(sorted_var_Sex_Cell)
#transpose back to rows subjects and columns genes for PCA and get rid of row_var
t_var_Sex_EV_Cell <- as.data.frame(t(sorted_var_Sex_EV_Cell[,-37]))
PCA_var_Sex_EV_Cell<-prcomp(t_var_Sex_EV_Cell, center = TRUE,scale. = TRUE)
t_var_Sex_EV <- as.data.frame(t(sorted_var_Sex_EV[,-19]))
PCA_var_Sex_EV<-prcomp(t_var_Sex_EV, center = TRUE,scale. = TRUE)
t_var_Sex_Cell <- as.data.frame(t(sorted_var_Sex_Cell[,-19]))
PCA_var_Sex_Cell<-prcomp(t_var_Sex_Cell, center = TRUE,scale. = TRUE)
View(PCA_var_Sex_Cell)

#samples grouped by EV and Cell
PCA_var_Sex_EV_Cell_components <- as.data.frame(PCA_var_Sex_EV_Cell$rotation)
View(PCA_var_Sex_EV_Cell_components)
write.csv(PCA_var_Sex_EV_Cell_components, "Sex_EV_Cell PCA contributors.csv")

PCA_var_Sex_EV_components <- as.data.frame(PCA_var_Sex_EV$rotation)
write.csv(PCA_var_Sex_EV_components, "Sex_EV PCA contributors.csv")

PCA_var_Sex_Cell_components <- as.data.frame(PCA_var_Sex_Cell$rotation)
write.csv(PCA_var_Sex_Cell_components, "Sex_Cell PCA contributors.csv")


dtp <- data.frame('EV_vs_Cell' = tSex_EV_Cell_norm_counts$EV_vs_Cell, 
                  'Sex' = tSex_EV_Cell_norm_counts$Sex,
                  'Pregnancy' = tSex_EV_Cell_norm_counts$Pregnancy,
                  'Treatment' = tSex_EV_Cell_norm_counts$Treatment,
                  PCA_var_Sex_EV_Cell$x[,1:2]) # the first two components are selected (NB: you can also select 3 for 3D plottings or 3+)
View(dtp)

dtp_EV <- data.frame( 
  'Sex' = tSex_EV_norm_counts$Sex,
  'Pregnancy' = tSex_EV_norm_counts$Pregnancy,
  'Treatment' = tSex_EV_norm_counts$Treatment,
  PCA_var_Sex_EV$x[,1:2]) # the first two components are selected (NB: you can also select 3 for 3D plottings or 3+)
View(dtp_EV)

dtp_Cell <- data.frame(
  'Sex' = tSex_Cell_norm_counts$Sex,
  'Pregnancy' = tSex_Cell_norm_counts$Pregnancy,
  'Treatment' = tSex_Cell_norm_counts$Treatment,
  PCA_var_Sex_Cell$x[,1:2]) # the first two components are selected (NB: you can also select 3 for 3D plottings or 3+)
View(dtp_Cell)

#by EV_vs_Cell
png("Sex_EV_Cell_PCA_top_500_variant_by_EV_vs_Cell.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp) + 
  geom_point(aes(x = PC1, y = PC2, col = Sex), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = EV_vs_Cell)) +
  scale_color_manual(values = c("#D55E00", "#000000", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00"))
dev.off()

autoplot(PCA_var_Sex_EV_Cell)
#autoplot shows the percentage of pc1 and pc2, which were 52.98 and 9.26 in this case.
#so graph with those in x and y labels 

png("Sex_EV_Cell_PCA_top_500_variant_by_EV_vs_Cell.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp) + 
  geom_point(aes(x = PC1, y = PC2, col = Sex), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = EV_vs_Cell)) +
  xlab("PC 1 (52.98%)") + 
  ylab("PC 2 (9.26%)") +
  scale_color_manual(values = c("#D55E00", "#000000", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00"))
dev.off()

#by Sex
png("Sex_EV_Cell_PCA_top_500_variant_by_Sex.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp) + 
  geom_point(aes(x = PC1, y = PC2, col = EV_vs_Cell), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = Sex)) +
  scale_color_manual(values = c("#D55E00", "#000000", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00"))
dev.off()

#EV samples only
#by Sex
autoplot(PCA_var_Sex_EV)
autoplot(PCA_var_Sex_EV, loadings=TRUE, loadings.colour='darkorchid4', loadings.label=TRUE, loadings.label.size=3)

png("Sex_EV_PCA_top_500_variant_by_Sex.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp_EV) + 
  geom_point(aes(x = PC1, y = PC2, col = Sex), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = Sex)) +
  xlab("PC 1 (39.69%)") + 
  ylab("PC 2 (16.65%)") +
  scale_color_manual(values = c("#D55E00", "#000000", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00"))
dev.off()


png("Sex_EV_PCA_top_500_variant_by_Pregnancy_ellipse.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp_EV) + 
  geom_point(aes(x = PC1, y = PC2, col = Pregnancy), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = Pregnancy)) +
  xlab("PC 1 (39.69%)") + 
  ylab("PC 2 (16.65%)") +
  scale_color_manual(values = c("#999999", "#009E73", "#F0E442", "#CC79A7", "#0072B2", "#E69F00"))
dev.off()

png("Sex_EV_PCA_top_500_variant_by_Pregnancy_ellipse.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp_EV) + 
  geom_point(aes(x = PC1, y = PC2, col = Pregnancy), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = Pregnancy)) +
  xlab("PC 1 (39.69%)") + 
  ylab("PC 2 (16.65%)") +
  scale_color_manual(values = c("#56B4E9", "#009E73", "#CC79A7"))
dev.off()

png("Sex_EV_PCA_top_500_variant_by_Sex_ellipse.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp_EV) + 
  geom_point(aes(x = PC1, y = PC2, col = Pregnancy), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = Sex)) +
  xlab("PC 1 (39.69%)") + 
  ylab("PC 2 (16.65%)") +
  scale_color_manual(values = c("#999999", "#009E73", "#F0E442", "#CC79A7", "#0072B2", "#E69F00"))
dev.off()

+
  scale_y_reverse()+ scale_x_reverse()

#Cell samples only
#by Sex
autoplot(PCA_var_Sex_Cell)
autoplot(PCA_var_Sex_Cell, loadings=TRUE, loadings.colour='darkorchid4', loadings.label=TRUE, loadings.label.size=3)

png("Sex_Cell_PCA_top_500_variant_by_Sex.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp_Cell) + 
  geom_point(aes(x = PC1, y = PC2, col = Sex), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = Sex)) +
  xlab("PC 1 (47.98%)") + 
  ylab("PC 2 (17.4%)") +
  scale_color_manual(values = c("#D55E00", "#000000", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00"))
dev.off()


png("Sex_Cell_PCA_top_500_variant_by_Sex_ellipse.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp_Cell) + 
  geom_point(aes(x = PC1, y = PC2, col = Pregnancy), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = Sex)) +
  xlab("PC 1 (47.98%)") + 
  ylab("PC 2 (17.4%)") +
  scale_color_manual(values = c("#D55E00", "#000000", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00"))
dev.off()

png("Sex_Cell_PCA_top_500_variant_by_Pregnancy_ellipse.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp_Cell) + 
  geom_point(aes(x = PC1, y = PC2, col = Pregnancy), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = Pregnancy)) +
  xlab("PC 1 (47.98%)") + 
  ylab("PC 2 (17.4%)") +
  scale_color_manual(values = c("#56B4E9", "#009E73", "#CC79A7"))
dev.off()

png("Sex_Cell_PCA_top_500_variant_by_Sex_ellipse.png", width = 2500, height = 2500, res = 300)
ggplot(data = dtp_Cell) + 
  geom_point(aes(x = PC1, y = PC2, col = Pregnancy), size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, col = Sex)) +
  xlab("PC 1 (47.98%)") + 
  ylab("PC 2 (17.4%)") +
  scale_color_manual(values = c("#999999", "#009E73", "#F0E442", "#CC79A7", "#0072B2", "#E69F00"))
dev.off()




####Volcano Plot of EV_vs_Cell_results
#https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

install.packages("ggplot2")
install.packages("ggrepel")

library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)


View(Sex_EV_Cell_results)
#The default cut-off for log2FC is >|1|; the default cut-off for P value is 10e-6. 
png("Sex_EV_Cell_Volcano_plot_default.png", width = 2500, height = 2500, res = 300)
EnhancedVolcano(Sex_EV_Cell_results,
                lab = rownames(Sex_EV_Cell_results),
                title = 'Female vs Male',
                x = 'log2(FC)',
                y = 'P-adj')
dev.off()

#custom cut-offs: FC_log2 0.5-1.5; default adj_pval of 10e-6
png("ANOVA_Sex_EV_Cell_Volcano_plot_zoom.png", width = 2500, height = 2500, res = 300)
EnhancedVolcano(Sex_EV_Cell_results,
                lab = rownames(Sex_EV_Cell_results),
                x = 'log2(FC)',
                y = 'P-adj',
                xlim = c(-2.5, 2.5),
                ylim = c(0, 15),
                title = 'Female vs Male',
                pCutoff = 10e-6,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30')
dev.off()

png("ANOVA_Sex_EV_Cell_Volcano_plot_zoom_10e-3.png", width = 2500, height = 2500, res = 300)
EnhancedVolcano(Sex_EV_Cell_results,
                lab = rownames(Sex_EV_Cell_results),
                x = 'log2(FC)',
                y = 'P-adj',
                title = 'Female vs Male',
                pCutoff = 10e-3,
                xlim = c(-3, 3),
                ylim = c(0, 15))
dev.off()

View(Sex_EV_results)

png("Sex_EV_Volcano_plot_default.png", width = 2500, height = 2500, res = 300)
EnhancedVolcano(Sex_EV_results,
                lab = rownames(Sex_EV_results),
                title = 'EV Female vs Male',
                x = 'log2(FC)',
                y = 'P-adj')
dev.off()

#let's zoom in, which will leave some genes out
#just put both plots together
png("Sex_EV_Volcano_plot_padj_10e-3_zoom.png", width = 2500, height = 2500, res = 300)
EnhancedVolcano(Sex_EV_results,
                lab = rownames(Sex_EV_results),
                title = 'EV Female vs Male',
                x = 'log2(FC)',
                y = 'P-adj',
                pCutoff = 10e-3,
                xlim = c(-3, 3),
                ylim = c(0, 15))
dev.off()

pdf("Sex_EV_Volcano_plot_default.pdf")
EnhancedVolcano(Sex_EV_results,
                lab = rownames(Sex_EV_results),
                title = 'EV Female vs Male',
                x = 'log2(FC)',
                y = 'P-adj',
                xlim = c(-4, 4))
dev.off()

#custom cut-offs: FC_log2 0.5-1.5; default adj_pval of 10e-6
pdf("ANOVA_Sex_EV_Volcano_plot_custom.pdf")
EnhancedVolcano(Sex_EV_results,
                lab = rownames(Sex_EV_results),
                x = 'log2(FC)',
                y = 'P-adj',
                xlim = c(-2, 2),
                title = 'EV Female vs Male',
                pCutoff = 10e-6,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30')
dev.off()


png("Sex_Cell_Volcano_plot_default.png", width = 2500, height = 2500, res = 300)
EnhancedVolcano(Sex_Cell_results,
                lab = rownames(Sex_Cell_results),
                title = 'Cell Female vs Male',
                x = 'log2(FC)',
                y = 'P-adj')
dev.off()

#let's zoom in, which will leave 4 genes out
#just put both plots together
png("Sex_Cell_Volcano_plot_zoom_10e-3.png", width = 2500, height = 2500, res = 300)
EnhancedVolcano(Sex_Cell_results,
                lab = rownames(Sex_Cell_results),
                title = 'Cell Female vs Male',
                x = 'log2(FC)',
                y = 'P-adj',
                pCutoff = 10e-3,
                xlim = c(-3, 3),
                ylim = c(0, 15))
dev.off()




#previously from ev vs cell.
View(EV_vs_Cell_results)
png("EV_vs_Cell_Volcano_plot_default.png", width = 2500, height = 2500, res = 300)
EnhancedVolcano(EV_vs_Cell_results,
                lab = rownames(EV_vs_Cell_results),
                title = 'EV vs. Cell',
                x = 'log2(FC)',
                y = 'P-adj',
                xlim = c(-6, 6))
dev.off()


#maxoverlapsConnectors = Inf is optional and too many texts
getwd()
setwd("C:/Users/daehy/Dropbox/Transcriptomic/R")


########

####let's do pathway analysis for Sex_EV_results_sig_genes and Sex_Cell_results_sig_genes, adjusted pvalue of 0.05####
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(DOSE)

library(ReactomePA)
View(Sex_EV_results_sig_genes)
Sex_EV_results_sig_genes$Gene_ID<-rownames(Sex_EV_results_sig_genes)
ID_location <- EV_vs_Cell_results_sig_genes$Gene_ID
eg_location <- bitr(ID_location, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg_location)
#  SYMBOL         ENTREZID
#1      Tmem132e    270893
#2 1500009L16Rik     69784
#3          Abi3     66610
#4  LOC118568032 118568032
#5 9030624G23Rik     66808
#6       Gm30068 102631829


#rename column in eg so matching title
View(eg_location)
colnames(eg_location)[1] <- "Gene_ID"

#merge based on Gene_ID for log2(FC)
Sex_EV_results_Entrez_log2_FC <- 
  merge(eg_location, Sex_EV_results_sig_genes[,c(2,7)], by = "Gene_ID")
View(Sex_EV_results_Entrez_log2_FC)
#get rid of gene id so i have entrez and fc for 2 columns only
Sex_EV_results_Entrez_log2_FC <- Sex_EV_results_Entrez_log2_FC[,-1]
#make it numeric like geneList
#https://bioinformatics.stackexchange.com/questions/13036/reactomepa-dataset
myGeneList <- Sex_EV_results_Entrez_log2_FC[,2]
names(myGeneList) <- as.character(Sex_EV_results_Entrez_log2_FC[,1])
myGeneList <- sort(myGeneList, decreasing = TRUE)
head(myGeneList)
View(myGeneList)

#start processing for Reactome
#limiting genes with log2 fold change that are greater than 0.
#If we use log2(fold change), fold changes lower than 1 (when B > A) become negative, 
#while those greater than 1 (A > B) become positive. 
#Now the values are symmetrical and it's easier to see fold changes in both directions on one plot.
#so, we will do anything above 0 as enriched in EVs since that is fold change of greater than 1. 
#use abs() if you dont care if the fc is negative or positive.

#data(geneList, package="DOSE")
#View(geneList)
class(myGeneList) #its numeric
de <- names(myGeneList)[myGeneList > 0]
head(de)

#The enrichPathway function allows users to select an appropriate background of genes as the baseline. 
#The gsePathway function supports GSEA to evaluate the enriched reactome pathways of high-throughput data. 
#These approaches can be used to verify interesting altered pathways and to identify unanticipated pathway associations.
x <- enrichPathway(gene=de, organism = "mouse", pvalueCutoff=0.05, pAdjustMethod = "BH", readable=T)
head(as.data.frame(x))
er <- as.data.frame(x)
#enriched term found. when i did it with pc1, i got no result, which is a valid result. 
#getting enriched terms are good too. 

getwd()
write.csv(er, "DESEQ2_EV_results_Female_enriched_foldchange2_Reactome_Results.csv")

png("DESEQ2_EV_results_Female_enriched_foldchange2_bar_Reactome.png", width = 4000, height = 2000, res = 300)
barplot(x, showCategory=10)
dev.off()

png("DESEQ2_EV_results_Female_enriched_foldchange2_dotplot_Reactome.png", width = 2000, height = 2000, res = 300)
dotplot(x, showCategory=10)
dev.off()

png("DESEQ2_EV_results_Female_enriched_foldchange2_dotplot_Reactome.png", width = 1500, height = 2500, res = 300)
dotplot(x, showCategory=15, font.size = 10)
dev.off()



##https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
##additional dot plot adjustments!!!!!!

x2<- pairwise_termsim(x) 
emapplot(x2, color="pvalue")

png("DESEQ2_EV_results_Female_enriched_foldchange2__enrichMap_Reactome.png", width = 2000, height = 2000, res = 300)
emapplot(x2)
dev.off()

png("DESEQ2_EV_results_Female_enriched_foldchange2_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         cex_label_gene = 0.8)
dev.off()


png("DESEQ2_EV_results_Female_enriched_foldchange2__cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         node_label="category")
dev.off()

#GSEA
#9.3 Reactome pathway gene set enrichment analysis; gsePathway()
#Gene Set Enrichment Analysis (GSEA) is a computational method 
#that determines whether an a priori defined set of genes shows statistically
#significant, concordant differences between two biological states (e.g. phenotypes).
y <- gsePathway(myGeneList, 
                organism = "mouse",
                pvalueCutoff=0.05,
                pAdjustMethod="BH", verbose=FALSE)
y2 <- pairwise_termsim(y) 

png("DESEQ2_EV_results_Female_enriched_foldchange2_GSEA_Reactome.png", width = 3000, height = 2000, res = 300)
emapplot(y2)
dev.off()

#gsePathway works, because it has genes enriched in EVs or cells in the "myGeneList". 
#whenever enrichment score and NES scores are negative, 
#that means they are genes that are at the bottom of the list or those with negative log2(FC) scores for my case.
#this means i have results for over-represented pathways for EV enriched genes and for cell enriched genes. 

head(as.data.frame(y2))
View(as.data.frame(y))
head(as.data.frame(y2))
er <- as.data.frame(y2)

getwd()
write.csv(er, "DESEQ2_EV_results_Female_enriched_foldchange2_GSEA_Reactome.csv")

#bar and dotplots are useless for looking at ev enriched genes unless i set the order by enrichment score,
#where positive is ev enriched.
png("DESEQ2_EV_results_Female_enriched_foldchange2_bar_GSEA_Reactome.png", width = 4000, height = 2000, res = 300)
barplot(y2, showCategory=10)
dev.off()


png("DESEQ2_EV_results_Female_enriched_foldchange2_dotplot_GSEA_Reactome.png", width = 1500, height = 2500, res = 300)
dotplot(y2, showCategory=15, font.size = 10)
dev.off()

#KEGG enrichment analysis
#visualize KEGG

#7.2 KEGG pathway over-representation analysis


View(de)
kk <- enrichKEGG(gene         = de,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
head(kk)

View(as.data.frame(kk))
er_kk<- as.data.frame(kk)
write.csv(er_kk_locations, "DESEQ2_EV_results_Female_enriched_KEGG_Results.csv")


##let's try EV male enriched:####
de_ev_male <- names(myGeneList)[myGeneList < 0]

head(de_ev_male)
View(de_ev_male)
#de_ev_male<-rev(de_ev_male) do not need to reverse row order; does not affect the pathways analysis

#The enrichPathway function allows users to select an appropriate background of genes as the baseline. 
#The gsePathway function supports GSEA to evaluate the enriched reactome pathways of high-throughput data. 
#These approaches can be used to verify interesting altered pathways and to identify unanticipated pathway associations.
x <- enrichPathway(gene=de_ev_male, organism = "mouse", pvalueCutoff=0.05, pAdjustMethod = "BH", readable=T)
head(as.data.frame(x))
er <- as.data.frame(x)
#enriched term found for cells. 

getwd()
write.csv(er, "DESEQ2_EV_results_male_enriched_foldchange2_Reactome_Results.csv")

png("DESEQ2_EV_results_male_enriched_foldchange2_barplot_Reactome.png", width = 4000, height = 2000, res = 300)
barplot(x, showCategory=10)
dev.off()

png("DESEQ2_EV_results_male_enriched_foldchange2_dotplot_Reactome.png", width = 2000, height = 2000, res = 300)
dotplot(x, showCategory=10)
dev.off()

png("DESEQ2_EV_results_male_enriched_foldchange2_dotplot_Reactome.png", width = 1500, height = 2500, res = 300)
dotplot(x, showCategory=15, font.size = 10)
dev.off()



##https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
##additional dot plot adjustments!!!!!!

x2 <- pairwise_termsim(x) 
emapplot(x2, color="pvalue")

png("DESEQ2_EV_results_male_enriched_foldchange2_enrichMap_Reactome.png", width = 2000, height = 2000, res = 300)
emapplot(x2)
dev.off()

png("DESEQ2_EV_results_male_enriched_foldchange2_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         cex_label_gene = 0.8)
dev.off()

png("DESEQ2_EV_results_male_enriched_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         node_label="category")
dev.off()

#GSEA


#KEGG enrichment analysis
#visualize KEGG

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("pathview")
#library("pathview")

#7.2 KEGG pathway over-representation analysis


View(de_ev_male)
kk <- enrichKEGG(gene         = de_ev_male,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
head(kk)

View(as.data.frame(kk))
er_kk<- as.data.frame(kk)
getwd()
write.csv(er_kk, "DESEQ2_EV_results_male_enriched_foldchange2_KEGG_Results.csv")

#browseKEGG(kk, 'mmu03040')
#################

#let's do Cell female enriched

View(Sex_Cell_results_sig_genes)
Sex_Cell_results_sig_genes$Gene_ID<-rownames(Sex_Cell_results_sig_genes)
ID_location <- Sex_Cell_results_sig_genes$Gene_ID
eg_location <- bitr(ID_location, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg_location)
#  SYMBOL         ENTREZID
#1      Tmem132e    270893
#2 1500009L16Rik     69784
#3          Abi3     66610
#4  LOC118568032 118568032
#5 9030624G23Rik     66808
#6       Gm30068 102631829


#rename column in eg so matching title
View(eg_location)
colnames(eg_location)[1] <- "Gene_ID"

#merge based on Gene_ID for log2(FC)
Sex_Cell_results_Entrez_log2_FC <- 
  merge(eg_location, Sex_Cell_results_sig_genes[,c(2,7)], by = "Gene_ID")
View(Sex_Cell_results_Entrez_log2_FC)
#get rid of gene id so i have entrez and fc for 2 columns only
Sex_Cell_results_Entrez_log2_FC <- Sex_Cell_results_Entrez_log2_FC[,-1]
#make it numeric like geneList
#https://bioinformatics.stackexchange.com/questions/13036/reactomepa-dataset
myGeneList <- Sex_Cell_results_Entrez_log2_FC[,2]
names(myGeneList) <- as.character(Sex_Cell_results_Entrez_log2_FC[,1])
myGeneList <- sort(myGeneList, decreasing = TRUE)
head(myGeneList)
View(myGeneList)

#start processing for Reactome
#limiting genes with log2 fold change that are greater than 0.
#If we use log2(fold change), fold changes lower than 1 (when B > A) become negative, 
#while those greater than 1 (A > B) become positive. 
#Now the values are symmetrical and it's easier to see fold changes in both directions on one plot.
#so, we will do anything above 0 as enriched in EVs since that is fold change of greater than 1. 
#use abs() if you dont care if the fc is negative or positive.

#data(geneList, package="DOSE")
#View(geneList)
class(myGeneList) #its numeric
de <- names(myGeneList)[myGeneList > 0]
head(de)

#The enrichPathway function allows users to select an appropriate background of genes as the baseline. 
#The gsePathway function supports GSEA to evaluate the enriched reactome pathways of high-throughput data. 
#These approaches can be used to verify interesting altered pathways and to identify unanticipated pathway associations.
x <- enrichPathway(gene=de, organism = "mouse", pvalueCutoff=0.05, pAdjustMethod = "BH", readable=T)
head(as.data.frame(x))
er <- as.data.frame(x)
View(er)
#enriched term found. when i did it with pc1, i got no result, which is a valid result. 
#getting enriched terms are good too. 

getwd()
write.csv(er, "DESEQ2_Cell_results_Female_enriched_foldchange2_Reactome_Results.csv")

png("DESEQ2_Cell_results_Female_enriched_foldchange2_bar_Reactome.png", width = 4000, height = 2000, res = 300)
barplot(x, showCategory=10)
dev.off()

png("DESEQ2_Cell_results_Female_enriched_foldchange2_dotplot_Reactome.png", width = 2000, height = 2000, res = 300)
dotplot(x, showCategory=10)
dev.off()

png("DESEQ2_Cell_results_Female_enriched_foldchange2__dotplot_Reactome.png", width = 1500, height = 2500, res = 300)
dotplot(x, showCategory=15, font.size = 10)
dev.off()



##https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
##additional dot plot adjustments!!!!!!

x2<- pairwise_termsim(x) 
emapplot(x2, color="pvalue")

png("DESEQ2_Cell_results_Female_enriched_foldchange2__enrichMap_Reactome.png", width = 2000, height = 2000, res = 300)
emapplot(x2)
dev.off()

png("DESEQ2_Cell_results_Female_enriched_foldchange2_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         cex_label_gene = 0.8)
dev.off()


png("DESEQ2_Cell_results_Female_enriched_foldchange2__cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         node_label="category")
dev.off()

#GSEA
#9.3 Reactome pathway gene set enrichment analysis; gsePathway()
#Gene Set Enrichment Analysis (GSEA) is a computational method 
#that determines whether an a priori defined set of genes shows statistically
#significant, concordant differences between two biological states (e.g. phenotypes).
y <- gsePathway(myGeneList, 
                organism = "mouse",
                pvalueCutoff=0.05,
                pAdjustMethod="BH", verbose=FALSE)
y2 <- pairwise_termsim(y) 

png("DESEQ2_Cell_results_Female_enriched_foldchange2_GSEA_Reactome.png", width = 3000, height = 2000, res = 300)
emapplot(y2)
dev.off()

#gsePathway works, because it has genes enriched in EVs or cells in the "myGeneList". 
#whenever enrichment score and NES scores are negative, 
#that means they are genes that are at the bottom of the list or those with negative log2(FC) scores for my case.
#this means i have results for over-represented pathways for EV enriched genes and for cell enriched genes. 

head(as.data.frame(y2))
View(as.data.frame(y))
head(as.data.frame(y2))
er <- as.data.frame(y2)

getwd()
write.csv(er, "DESEQ2_Cell_results_Female_enriched_foldchange2_GSEA_Reactome.csv")

#bar and dotplots are useless for looking at ev enriched genes unless i set the order by enrichment score,
#where positive is ev enriched.
png("DESEQ2_Cell_results_Female_enriched_foldchange2_bar_GSEA_Reactome.png", width = 4000, height = 2000, res = 300)
barplot(y2, showCategory=10)
dev.off()


png("DESEQ2_Cell_results_Female_enriched_foldchange2_dotplot_GSEA_Reactome.png", width = 1500, height = 2500, res = 300)
dotplot(y2, showCategory=15, font.size = 10)
dev.off()

#KEGG enrichment analysis
#visualize KEGG

#7.2 KEGG pathway over-representation analysis


View(de)
kk <- enrichKEGG(gene         = de,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
head(kk)

View(as.data.frame(kk))
er_kk<- as.data.frame(kk)
write.csv(er_kk_locations, "DESEQ2_Cell_results_Female_enriched_KEGG_Results.csv")


##let's try Cell male enriched:####
de_cell_male <- names(myGeneList)[myGeneList < 0]

head(de_cell_male)
View(de_cell_male)
#de_cell_male<-rev(de_cell_male) do not need to reverse row order; does not affect the pathways analysis

#The enrichPathway function allows users to select an appropriate background of genes as the baseline. 
#The gsePathway function supports GSEA to evaluate the enriched reactome pathways of high-throughput data. 
#These approaches can be used to verify interesting altered pathways and to identify unanticipated pathway associations.
x <- enrichPathway(gene=de_cell_male, organism = "mouse", pvalueCutoff=0.05, pAdjustMethod = "BH", readable=T)
head(as.data.frame(x))
er <- as.data.frame(x)
#enriched term found for cells. 

getwd()
write.csv(er, "DESEQ2_Cell_results_male_enriched_foldchange2_Reactome_Results.csv")

png("DESEQ2_Cell_results_male_enriched_foldchange2_barplot_Reactome.png", width = 4000, height = 2000, res = 300)
barplot(x, showCategory=10)
dev.off()

png("DESEQ2_Cell_results_male_enriched_foldchange2_dotplot_Reactome.png", width = 2000, height = 2000, res = 300)
dotplot(x, showCategory=10)
dev.off()

png("DESEQ2_Cell_results_male_enriched_foldchange2__dotplot_Reactome.png", width = 1500, height = 2500, res = 300)
dotplot(x, showCategory=15, font.size = 10)
dev.off()



##https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
##additional dot plot adjustments!!!!!!

x2 <- pairwise_termsim(x) 
emapplot(x2, color="pvalue")

png("DESEQ2_Cell_results_male_enriched_foldchange2_enrichMap_Reactome.png", width = 2000, height = 2000, res = 300)
emapplot(x2)
dev.off()

png("DESEQ2_Cell_results_male_enriched_foldchange2_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         cex_label_gene = 0.8)
dev.off()

png("DESEQ2_Cell_results_male_enriched_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x2, categorySize="pvalue", foldChange=myGeneList,
         cex_label_category = 0.8, 
         node_label="category")
dev.off()

#GSEA


#KEGG enrichment analysis
#visualize KEGG

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("pathview")
#library("pathview")

#7.2 KEGG pathway over-representation analysis


View(de_cell_male)
kk <- enrichKEGG(gene         = de_cell_male,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
head(kk)

View(as.data.frame(kk))
er_kk<- as.data.frame(kk)
getwd()
write.csv(er_kk, "DESEQ2_Cell_results_male_enriched_foldchange2_KEGG_Results.csv")

###################




#####etoh effsize and value####
#let's do paired t test and effect size for ethanol treatment
library(effsize)

View(EV_vs_Cell_norm_counts_Gene_ID)
View(t_EV_vs_Cell_norm_counts_Gene_ID)
View(n_ev22_cell22_NA_merge_then_normalize_byEVproteins_1)
View(ntall_n_ev22_cell22_NA_merge_then_normalize_byEVproteins_1)

rna_expression_ev_cell<-t_EV_vs_Cell_norm_counts_Gene_ID
#change columns from character to numeric
# Convert all variable types to numeric
rna_expression_ev_cell <- as.data.frame(apply(rna_expression_ev_cell, 2, as.numeric))  
# Print classes of all columns
#sapply(rna_expression_ev_cell, class)  
#i lose row names so let's put that back in
rownames(rna_expression_ev_cell)<-rownames(t_EV_vs_Cell_norm_counts_Gene_ID)

#rows are proteins and columns are samples. 

View(rna_expression_ev_cell)
tp_table_rna_ev_ctrl_vs_120<-data.frame(RNA = rep(NA,ncol(rna_expression_ev_cell)), 
                                        pvalue = rep(NA,ncol(rna_expression_ev_cell)),
                                        effsize=rep(NA,ncol(rna_expression_ev_cell)),
                                        mean_control=rep(NA,ncol(rna_expression_ev_cell)),
                                        mean_ethanol=rep(NA,ncol(rna_expression_ev_cell)))
View(tp_table_rna_ev_ctrl_vs_120)




#rna starts from column 1 in this case, hence i in 1:
for (i in 1:ncol(rna_expression_ev_cell))
{evttest<-t.test(rna_expression_ev_cell[c(1,4,7,10,13,16),i],
                 rna_expression_ev_cell[c(2,5,8,11,14,17),i],
                 paired=TRUE)
tp_table_rna_ev_ctrl_vs_120[i,2]<-evttest$p.value
tp_table_rna_ev_ctrl_vs_120[i,1]<-colnames(rna_expression_ev_cell)[i]
ef_size<-cohen.d(rna_expression_ev_cell[c(2,5,8,11,14,17),i],
                 rna_expression_ev_cell[c(1,4,7,10,13,16),i],
                 pooled=TRUE,paired=TRUE,na.rm=TRUE, hedges.correction=TRUE)
tp_table_rna_ev_ctrl_vs_120[i,3]<-ef_size$estimate
tp_table_rna_ev_ctrl_vs_120[i,5]<-mean(rna_expression_ev_cell[c(2,5,8,11,14,17),i])
tp_table_rna_ev_ctrl_vs_120[i,4]<-mean(rna_expression_ev_cell[c(1,4,7,10,13,16),i])}
dim(tp_table_rna_ev_ctrl_vs_120)
View(tp_table_rna_ev_ctrl_vs_120)
summary(tp_table_rna_ev_ctrl_vs_120$mean_control)
summary(tp_table_rna_ev_ctrl_vs_120$mean_ethanol)

#now 0vs320
tp_table_rna_ev_ctrl_vs_320<-data.frame(RNA = rep(NA,ncol(rna_expression_ev_cell)), 
                                        pvalue = rep(NA,ncol(rna_expression_ev_cell)),
                                        effsize=rep(NA,ncol(rna_expression_ev_cell)),
                                        mean_control=rep(NA,ncol(rna_expression_ev_cell)),
                                        mean_ethanol=rep(NA,ncol(rna_expression_ev_cell)))

#rna starts from column 1 in this case, hence i in 1:
for (i in 1:ncol(rna_expression_ev_cell))
{evttest<-t.test(rna_expression_ev_cell[c(1,4,7,10,13,16),i],
                 rna_expression_ev_cell[c(3,6,9,12,15,18),i],
                 paired=TRUE)
tp_table_rna_ev_ctrl_vs_320[i,2]<-evttest$p.value
tp_table_rna_ev_ctrl_vs_320[i,1]<-colnames(rna_expression_ev_cell)[i]
ef_size<-cohen.d(rna_expression_ev_cell[c(3,6,9,12,15,18),i],
                 rna_expression_ev_cell[c(1,4,7,10,13,16),i],
                 pooled=TRUE,paired=TRUE,na.rm=TRUE, hedges.correction=TRUE)
tp_table_rna_ev_ctrl_vs_320[i,3]<-ef_size$estimate
tp_table_rna_ev_ctrl_vs_320[i,5]<-mean(rna_expression_ev_cell[c(3,6,9,12,15,18),i])
tp_table_rna_ev_ctrl_vs_320[i,4]<-mean(rna_expression_ev_cell[c(1,4,7,10,13,16),i])}
dim(tp_table_rna_ev_ctrl_vs_320)
View(tp_table_rna_ev_ctrl_vs_320)
summary(tp_table_rna_ev_ctrl_vs_320$mean_control)
summary(tp_table_rna_ev_ctrl_vs_320$mean_ethanol)

#let's see how i can merge the data, 
#so I can have 0vs320 pvalue and effect size next to 0vs120 proteins with significant pvalue, and vice versa.
#better yet, lets merge the two d.f. together.
#to do that, let's rename some columns to specify which treatment it comes from. 
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))

library(data.table)

setnames(tp_table_rna_ev_ctrl_vs_120, old=c("pvalue","effsize","mean_control","mean_ethanol"), new=c("pvalue_0vs120", "effsize_0vs120","mean_control","mean_ethanol_120"))
setnames(tp_table_rna_ev_ctrl_vs_320, old=c("pvalue","effsize","mean_control","mean_ethanol"), new=c("pvalue_0vs320", "effsize_0vs320","mean_control","mean_ethanol_320"))
#for merging into a new d.f., let's only include columns I am interested in merging into. 
rna_ev_t_test_and_effect_size_all<-merge(tp_table_rna_ev_ctrl_vs_120,tp_table_rna_ev_ctrl_vs_320,by="Protein")

View(rna_ev_t_test_and_effect_size_all)

#reorder columns
#ev22_t_test_and_effect_size_all<-ev22_t_test_and_effect_size_all[,c(2,3,1,6,4,5,7,8,9,10)]
#View(ev22_t_test_and_effect_size_all)
#it worked but there are empty rows. remove that
#ev22_t_test_and_effect_size_all<-ev22_t_test_and_effect_size_all[1:2500,]


#let's add expression level to the list (detailed info of tSum below). 
#create a dataframe with a new row that is "Protein_Sum_ev22"
View(rna_expression_ev_cell)

rna_ev_Sum<-data.frame(matrix(NA, nrow = 18, ncol = 1))

#now, find the sum of each row
#practice on a subset; 1:18 rows are ev's. start from column 1 which is the first protein
View(rna_expression_ev_cell)
#rowSums(rna_expression_ev_cell[1:18,1:2], na.rm = TRUE)
colSums(rna_expression_ev_cell[1:18,1:2], na.rm = TRUE)

#worked. so now, lets find row sums then transfer that into the dataset "Sample_Sum_ev3" that i created
#Sample_Sum_ev5<-rowSums(rna_expression_ev_cell[1:18,6:2505], na.rm = TRUE)
rna_ev_Sum<-colSums(rna_expression_ev_cell[1:18,1:ncol(rna_expression_ev_cell)], na.rm = TRUE)
View(rna_ev_Sum)
median(rna_ev_Sum)
#scientific(median(rna_ev_Sum))
formatC(median(rna_ev_Sum),format = "e", digits = 2)
#49.00342
#"4.90e+01"


View(rna_ev_t_test_and_effect_size_all)
colnames(rna_ev_t_test_and_effect_size_all)[colnames(rna_ev_t_test_and_effect_size_all)=="Protein"] <- "Gene_ID"
rownames(rna_ev_t_test_and_effect_size_all)<-rna_ev_t_test_and_effect_size_all$Gene_ID
#merge
rna_ev_t_test_and_effect_size_all_Sum<-merge(rna_ev_t_test_and_effect_size_all, rna_ev_Sum, by="row.names", all.x = T)
View(rna_ev_t_test_and_effect_size_all_Sum)


#there is mean_ctrl.x and mean_ctrl.y with same values so delete one
rna_ev_t_test_and_effect_size_all_Sum<-rna_ev_t_test_and_effect_size_all_Sum[,-c(1,9)]


# Rename a column in R
colnames(rna_ev_t_test_and_effect_size_all_Sum)[colnames(rna_ev_t_test_and_effect_size_all_Sum)=="mean_control.x"] <- "mean_ctrl"
colnames(rna_ev_t_test_and_effect_size_all_Sum)[colnames(rna_ev_t_test_and_effect_size_all_Sum)=="y"] <- "Total_Normalized_rna_Expression_Level_EV"

View(rna_ev_t_test_and_effect_size_all_Sum)

#remove duplicates
#it worked. but there are duplicate rows, so let's remove that with the function duplicated() and unique() to extract only unique rows
library(tidyverse)
#If you want to remove duplicated elements, use !duplicated(), where ! is a logical negation, telling we don't want duplicate rows.
#Remove duplicates based on Accession Number column, since R doesnt recognize space, rename Accession Number to Accession_Number
#colnames(rna_ev_t_test_and_effect_size_all_Sum)[3]<-"Accession_Number"
rna_ev_t_test_and_effect_size_all_Sum<-rna_ev_t_test_and_effect_size_all_Sum %>% distinct(Gene_ID, .keep_all = TRUE)
View(rna_ev_t_test_and_effect_size_all_Sum)

write.csv(rna_ev_t_test_and_effect_size_all_Sum,file="rna_ev_t_test_and_effect_size_all_Sum.csv")

###########################################################




#Let's use DESeq2 in EV samples and in Cell samples,
#with the combined criteria of raw p-value < 0.05; effsize > 0.4 with nonzero-containing 95% CI
#for etoh treatments (0vs120 and 0vs320) in EV samples and in Cell samples. 

##let's try 0vs120 then 0vs320 by ev samples and cell samples separate. 
#then maybe even more separate by female samples and male samples separate. 



#####etoh effsize and value####
#let's do paired t test and effect size for all genes with nonzero base mean
library(effsize)
view(t_rna_EV_Cell_norm_counts_nonzero_base_mean)
View(t_EV_vs_Cell_norm_counts_Gene_ID)
View(n_ev22_cell22_NA_merge_then_normalize_byEVproteins_1)
View(ntall_n_ev22_cell22_NA_merge_then_normalize_byEVproteins_1)

rna_expression_ev_cell<-t_EV_vs_Cell_norm_counts_Gene_ID
#change columns from character to numeric
# Convert all variable types to numeric
rna_expression_ev_cell <- as.data.frame(apply(rna_expression_ev_cell, 2, as.numeric))  
# Print classes of all columns
#sapply(rna_expression_ev_cell, class)  
#i lose row names so let's put that back in
rownames(rna_expression_ev_cell)<-rownames(t_EV_vs_Cell_norm_counts_Gene_ID)

#rows are rnas and columns are samples. 
#The function computes the value of Cohen's d statistics (Cohen 1988). 
#If required (hedges.correction==TRUE) the Hedges g statistics is computed instead (Hedges and Holkin, 1985),
#which is what I did.


View(rna_expression_ev_cell)
tp_table_rna_ev_ctrl_vs_120<-data.frame(Gene_ID = rep(NA,ncol(rna_expression_ev_cell)), 
                                        pvalue = rep(NA,ncol(rna_expression_ev_cell)),
                                        effsize=rep(NA,ncol(rna_expression_ev_cell)),
                                        effsize_lower_CI=rep(NA,ncol(rna_expression_ev_cell)),
                                        effsize_upper_CI=rep(NA,ncol(rna_expression_ev_cell)),
                                        mean_control=rep(NA,ncol(rna_expression_ev_cell)),
                                        mean_ethanol=rep(NA,ncol(rna_expression_ev_cell)))
View(tp_table_rna_ev_ctrl_vs_120)




#gene starts from column 1 in this case, hence i in 1:
for (i in 1:ncol(rna_expression_ev_cell))
{evttest<-t.test(rna_expression_ev_cell[c(1,4,7,10,13,16),i],
                 rna_expression_ev_cell[c(2,5,8,11,14,17),i],
                 paired=TRUE)
tp_table_rna_ev_ctrl_vs_120[i,2]<-evttest$p.value
tp_table_rna_ev_ctrl_vs_120[i,1]<-colnames(rna_expression_ev_cell)[i]
ef_size<-cohen.d(rna_expression_ev_cell[c(2,5,8,11,14,17),i],
                 rna_expression_ev_cell[c(1,4,7,10,13,16),i],
                 pooled=TRUE,paired=TRUE,na.rm=TRUE, hedges.correction=TRUE)
tp_table_rna_ev_ctrl_vs_120[i,3]<-ef_size$estimate
tp_table_rna_ev_ctrl_vs_120[i,4]<-ef_size$conf.int[[1]]
tp_table_rna_ev_ctrl_vs_120[i,5]<-ef_size$conf.int[[2]]
tp_table_rna_ev_ctrl_vs_120[i,7]<-mean(rna_expression_ev_cell[c(2,5,8,11,14,17),i])
tp_table_rna_ev_ctrl_vs_120[i,6]<-mean(rna_expression_ev_cell[c(1,4,7,10,13,16),i])}
dim(tp_table_rna_ev_ctrl_vs_120)
View(tp_table_rna_ev_ctrl_vs_120)
summary(tp_table_rna_ev_ctrl_vs_120$mean_control)
summary(tp_table_rna_ev_ctrl_vs_120$mean_ethanol)

#now 0vs320
tp_table_rna_ev_ctrl_vs_320<-data.frame(Gene_ID = rep(NA,ncol(rna_expression_ev_cell)), 
                                        pvalue = rep(NA,ncol(rna_expression_ev_cell)),
                                        effsize=rep(NA,ncol(rna_expression_ev_cell)),
                                        effsize_lower_CI=rep(NA,ncol(rna_expression_ev_cell)),
                                        effsize_upper_CI=rep(NA,ncol(rna_expression_ev_cell)),
                                        mean_control=rep(NA,ncol(rna_expression_ev_cell)),
                                        mean_ethanol=rep(NA,ncol(rna_expression_ev_cell)))

#protein starts from column 1 in this case, hence i in 1:
for (i in 1:ncol(rna_expression_ev_cell))
{evttest<-t.test(rna_expression_ev_cell[c(1,4,7,10,13,16),i],
                 rna_expression_ev_cell[c(3,6,9,12,15,18),i],
                 paired=TRUE)
tp_table_rna_ev_ctrl_vs_320[i,2]<-evttest$p.value
tp_table_rna_ev_ctrl_vs_320[i,1]<-colnames(rna_expression_ev_cell)[i]
ef_size<-cohen.d(rna_expression_ev_cell[c(3,6,9,12,15,18),i],
                 rna_expression_ev_cell[c(1,4,7,10,13,16),i],
                 pooled=TRUE,paired=TRUE,na.rm=TRUE, hedges.correction=TRUE)
tp_table_rna_ev_ctrl_vs_320[i,3]<-ef_size$estimate
tp_table_rna_ev_ctrl_vs_320[i,4]<-ef_size$conf.int[[1]]
tp_table_rna_ev_ctrl_vs_320[i,5]<-ef_size$conf.int[[2]]
tp_table_rna_ev_ctrl_vs_320[i,7]<-mean(rna_expression_ev_cell[c(3,6,9,12,15,18),i])
tp_table_rna_ev_ctrl_vs_320[i,6]<-mean(rna_expression_ev_cell[c(1,4,7,10,13,16),i])}
dim(tp_table_rna_ev_ctrl_vs_320)
View(tp_table_rna_ev_ctrl_vs_320)
summary(tp_table_rna_ev_ctrl_vs_320$mean_control)
summary(tp_table_rna_ev_ctrl_vs_320$mean_ethanol)

#let's see how i can merge the data, 
#so I can have 0vs320 pvalue and effect size next to 0vs120 proteins with significant pvalue, and vice versa.
#better yet, lets merge the two d.f. together.
#to do that, let's rename some columns to specify which treatment it comes from. 
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))

library(data.table)

setnames(tp_table_rna_ev_ctrl_vs_120, old=c("pvalue","effsize","effsize_lower_CI","effsize_upper_CI","mean_control","mean_ethanol"), 
         new=c("pvalue_0vs120", "effsize_0vs120","effsize_0vs120_lower_CI","effsize_0vs120_upper_CI","mean_control","mean_ethanol_120"))
setnames(tp_table_rna_ev_ctrl_vs_320, old=c("pvalue","effsize","effsize_lower_CI","effsize_upper_CI","mean_control","mean_ethanol"), 
         new=c("pvalue_0vs320", "effsize_0vs320","effsize_0vs320_lower_CI","effsize_0vs320_upper_CI","mean_control","mean_ethanol_320"))
#for merging into a new d.f., let's only include columns I am interested in merging into. 
rna_ev_t_test_and_effect_size_all<-merge(tp_table_rna_ev_ctrl_vs_120,tp_table_rna_ev_ctrl_vs_320,by="Gene_ID")

View(rna_ev_t_test_and_effect_size_all)



#let's add expression level to the list (detailed info of tSum below). 
#create a dataframe with a new row that is "rna_ev_Sum"

rna_ev_Sum<-data.frame(matrix(NA, nrow = 18, ncol = 1))

#now, find the sum of each row
#practice on a subset; 1:18 rows are ev's. start from column 1 which is the first rna
View(rna_expression_ev_cell)
#rowSums(rna_expression_ev_cell[1:18,1:2], na.rm = TRUE)
colSums(rna_expression_ev_cell[1:18,1:2], na.rm = TRUE)

#worked. so now, lets find row sums then transfer that into the dataset 
rna_ev_Sum<-colSums(rna_expression_ev_cell[1:18,1:ncol(rna_expression_ev_cell)], na.rm = TRUE)
View(as.data.frame(rna_ev_Sum))
median(rna_ev_Sum)
#scientific(median(rna_ev_Sum))
formatC(median(rna_ev_Sum),format = "e", digits = 2)
#49.00342
#"4.90e+01"

View(rna_ev_t_test_and_effect_size_all)
rownames(rna_ev_t_test_and_effect_size_all)<-rna_ev_t_test_and_effect_size_all$Gene_ID
#merge
rna_ev_t_test_and_effect_size_all_Sum<-merge(rna_ev_t_test_and_effect_size_all, rna_ev_Sum, by="row.names", all.x = T)
View(rna_ev_t_test_and_effect_size_all_Sum)


#there is mean_ctrl.x and mean_ctrl.y with same values so delete one
rna_ev_t_test_and_effect_size_all_Sum<-rna_ev_t_test_and_effect_size_all_Sum[,-c(1,13)]


# Rename a column in R
colnames(rna_ev_t_test_and_effect_size_all_Sum)[colnames(rna_ev_t_test_and_effect_size_all_Sum)=="mean_control.x"] <- "mean_ctrl"
colnames(rna_ev_t_test_and_effect_size_all_Sum)[colnames(rna_ev_t_test_and_effect_size_all_Sum)=="y"] <- "Total_Normalized_rna_Expression_Level_EV"

View(rna_ev_t_test_and_effect_size_all_Sum)

#remove duplicates
#it worked. but there are duplicate rows, so let's remove that with the function duplicated() and unique() to extract only unique rows
library(tidyverse)
#If you want to remove duplicated elements, use !duplicated(), where ! is a logical negation, telling we don't want duplicate rows.
#Remove duplicates based on Accession Number column, since R doesnt recognize space, rename Accession Number to Accession_Number
#colnames(rna_ev_t_test_and_effect_size_all_Sum)[3]<-"Accession_Number"
rna_ev_t_test_and_effect_size_all_Sum<-rna_ev_t_test_and_effect_size_all_Sum %>% distinct(Gene_ID, .keep_all = TRUE)
View(rna_ev_t_test_and_effect_size_all_Sum)

write.csv(rna_ev_t_test_and_effect_size_all_Sum,file="rna_ev_t_test_and_effect_size_all_Sum.csv")



####now let's do cell's paired t test and effect size####
library(effsize)



tp_table_rna_Cell_ctrl_vs_120<-data.frame(Gene_ID = rep(NA,ncol(rna_expression_ev_cell)), 
                                          pvalue = rep(NA,ncol(rna_expression_ev_cell)),
                                          effsize=rep(NA,ncol(rna_expression_ev_cell)),
                                          effsize_lower_CI=rep(NA,ncol(rna_expression_ev_cell)),
                                          effsize_upper_CI=rep(NA,ncol(rna_expression_ev_cell)),
                                          mean_control=rep(NA,ncol(rna_expression_ev_cell)),
                                          mean_ethanol=rep(NA,ncol(rna_expression_ev_cell)))
View(tp_table_rna_Cell_ctrl_vs_120)




#protein starts from column 1 in this case, hence i in 1:
for (i in 1:ncol(rna_expression_ev_cell))
{evttest<-t.test(rna_expression_ev_cell[c(19,22,25,28,31,34),i],
                 rna_expression_ev_cell[c(20,23,26,29,32,35),i],
                 paired=TRUE)
tp_table_rna_Cell_ctrl_vs_120[i,2]<-evttest$p.value
tp_table_rna_Cell_ctrl_vs_120[i,1]<-colnames(rna_expression_ev_cell)[i]
ef_size<-cohen.d(rna_expression_ev_cell[c(20,23,26,29,32,35),i],
                 rna_expression_ev_cell[c(19,22,25,28,31,34),i],
                 pooled=TRUE,paired=TRUE,na.rm=TRUE, hedges.correction=TRUE)
tp_table_rna_Cell_ctrl_vs_120[i,3]<-ef_size$estimate
tp_table_rna_Cell_ctrl_vs_120[i,4]<-ef_size$conf.int[[1]]
tp_table_rna_Cell_ctrl_vs_120[i,5]<-ef_size$conf.int[[2]]
tp_table_rna_Cell_ctrl_vs_120[i,7]<-mean(rna_expression_ev_cell[c(20,23,26,29,32,35),i])
tp_table_rna_Cell_ctrl_vs_120[i,6]<-mean(rna_expression_ev_cell[c(19,22,25,28,31,34),i])}
dim(tp_table_rna_Cell_ctrl_vs_120)
View(tp_table_rna_Cell_ctrl_vs_120)
summary(tp_table_rna_Cell_ctrl_vs_120$mean_control)
summary(tp_table_rna_Cell_ctrl_vs_120$mean_ethanol)



#now 0vs320
tp_table_rna_Cell_ctrl_vs_320<-data.frame(Gene_ID = rep(NA,ncol(rna_expression_ev_cell)), 
                                          pvalue = rep(NA,ncol(rna_expression_ev_cell)),
                                          effsize=rep(NA,ncol(rna_expression_ev_cell)),
                                          effsize_lower_CI=rep(NA,ncol(rna_expression_ev_cell)),
                                          effsize_upper_CI=rep(NA,ncol(rna_expression_ev_cell)),
                                          mean_control=rep(NA,ncol(rna_expression_ev_cell)),
                                          mean_ethanol=rep(NA,ncol(rna_expression_ev_cell)))

#protein starts from column 1 in this case, hence i in 1:
for (i in 1:ncol(rna_expression_ev_cell))
{evttest<-t.test(rna_expression_ev_cell[c(19,22,25,28,31,34),i],
                 rna_expression_ev_cell[c(21,24,27,30,33,36),i],
                 paired=TRUE)
tp_table_rna_Cell_ctrl_vs_320[i,2]<-evttest$p.value
tp_table_rna_Cell_ctrl_vs_320[i,1]<-colnames(rna_expression_ev_cell)[i]
ef_size<-cohen.d(rna_expression_ev_cell[c(21,24,27,30,33,36),i],
                 rna_expression_ev_cell[c(19,22,25,28,31,34),i],
                 pooled=TRUE,paired=TRUE,na.rm=TRUE, hedges.correction=TRUE)
tp_table_rna_Cell_ctrl_vs_320[i,3]<-ef_size$estimate
tp_table_rna_Cell_ctrl_vs_320[i,4]<-ef_size$conf.int[[1]]
tp_table_rna_Cell_ctrl_vs_320[i,5]<-ef_size$conf.int[[2]]
tp_table_rna_Cell_ctrl_vs_320[i,7]<-mean(rna_expression_ev_cell[c(21,24,27,30,33,36),i])
tp_table_rna_Cell_ctrl_vs_320[i,6]<-mean(rna_expression_ev_cell[c(19,22,25,28,31,34),i])}
dim(tp_table_rna_Cell_ctrl_vs_320)
View(tp_table_rna_Cell_ctrl_vs_320)
summary(tp_table_rna_Cell_ctrl_vs_320$mean_control)
summary(tp_table_rna_Cell_ctrl_vs_320$mean_ethanol)


#let's see how i can merge the data, 
#so I can have 0vs320 pvalue and effect size next to 0vs120 proteins with significant pvalue, and vice versa.
#better yet, lets merge the two d.f. together.
#to do that, let's rename some columns to specify which treatment it comes from. 
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))

library(data.table)

setnames(tp_table_rna_Cell_ctrl_vs_120, old=c("pvalue","effsize","effsize_lower_CI","effsize_upper_CI","mean_control","mean_ethanol"), 
         new=c("pvalue_0vs120", "effsize_0vs120","effsize_0vs120_lower_CI","effsize_0vs120_upper_CI","mean_control","mean_ethanol_120"))
setnames(tp_table_rna_Cell_ctrl_vs_320, old=c("pvalue","effsize","effsize_lower_CI","effsize_upper_CI","mean_control","mean_ethanol"), 
         new=c("pvalue_0vs320", "effsize_0vs320","effsize_0vs320_lower_CI","effsize_0vs320_upper_CI","mean_control","mean_ethanol_320"))
#for merging into a new d.f., let's only include columns I am interested in merging into. 
rna_Cell_t_test_and_effect_size_all<-merge(tp_table_rna_Cell_ctrl_vs_120,tp_table_rna_Cell_ctrl_vs_320,by="Gene_ID")

View(rna_Cell_t_test_and_effect_size_all)

#reorder columns
#ev22_t_test_and_effect_size_all<-ev22_t_test_and_effect_size_all[,c(2,3,1,6,4,5,7,8,9,10)]
#View(ev22_t_test_and_effect_size_all)
#it worked but there are empty rows. remove that
#ev22_t_test_and_effect_size_all<-ev22_t_test_and_effect_size_all[1:2500,]


#let's add expression level to the list (detailed info of tSum below). 
#create a dataframe with a new row that is "Protein_Sum_ev22"


rna_Cell_Sum<-data.frame(matrix(NA, nrow = 18, ncol = 1))

#now, find the sum of each row
#practice on a subset; 1:18 rows are ev's. start from column 1 which is the first protein
View(rna_expression_ev_cell)
#rowSums(rna_expression_ev_cell[1:18,1:2], na.rm = TRUE)
colSums(rna_expression_ev_cell[1:18,1:2], na.rm = TRUE)

#worked. so now, lets find row sums then transfer that into the dataset "Sample_Sum_ev3" that i created
#Sample_Sum_ev5<-rowSums(rna_expression_ev_cell[1:18,6:2505], na.rm = TRUE)
rna_Cell_Sum<-colSums(rna_expression_ev_cell[19:36,1:ncol(rna_expression_ev_cell)], na.rm = TRUE)
View(as.data.frame(rna_Cell_Sum))
median(rna_Cell_Sum)
#scientific(median(rna_ev_Sum))
formatC(median(rna_Cell_Sum),format = "e", digits = 2)
#21.3168
#[1] "2.13e+01"
View(rna_Cell_Sum)
View(rna_Cell_t_test_and_effect_size_all)
View(as.data.frame(rna_Cell_Sum))
rownames(rna_Cell_t_test_and_effect_size_all)<-rna_Cell_t_test_and_effect_size_all$Gene_ID
#merge
rna_Cell_t_test_and_effect_size_all_Sum<-merge(rna_Cell_t_test_and_effect_size_all, rna_Cell_Sum, by="row.names", all.x = T)
View(rna_Cell_t_test_and_effect_size_all_Sum)


#there is mean_ctrl.x and mean_ctrl.y with same values so delete one
rna_Cell_t_test_and_effect_size_all_Sum<-rna_Cell_t_test_and_effect_size_all_Sum[,-c(1,13)]


#should merge with gene names and protein names
# Rename a column in R
colnames(rna_Cell_t_test_and_effect_size_all_Sum)[colnames(rna_Cell_t_test_and_effect_size_all_Sum)=="mean_control.x"] <- "mean_ctrl"
colnames(rna_Cell_t_test_and_effect_size_all_Sum)[colnames(rna_Cell_t_test_and_effect_size_all_Sum)=="y"] <- "Total_Normalized_rna_Expression_Level_Cell"

View(rna_Cell_t_test_and_effect_size_all_Sum)

#remove duplicates
#it worked. but there are duplicate rows, so let's remove that with the function duplicated() and unique() to extract only unique rows
library(tidyverse)
#If you want to remove duplicated elements, use !duplicated(), where ! is a logical negation, telling we don't want duplicate rows.
#Remove duplicates based on Accession Number column, since R doesnt recognize space, rename Accession Number to Accession_Number
#colnames(rna_ev_t_test_and_effect_size_all_Sum)[3]<-"Accession_Number"
rna_Cell_t_test_and_effect_size_all_Sum<-rna_Cell_t_test_and_effect_size_all_Sum %>% distinct(Gene_ID, .keep_all = TRUE)
View(rna_Cell_t_test_and_effect_size_all_Sum)

write.csv(rna_Cell_t_test_and_effect_size_all_Sum,file="rna_Cell_t_test_and_effect_size_all_Sum.csv")


#now, let's create a d.f. with only nonzero CI proteins.
#we can use subset() with multiple conditions, with & for "and", and | for "or".
View(rna_ev_t_test_and_effect_size_all_Sum)
#ev first#
ev_a120<-rna_ev_t_test_and_effect_size_all_Sum
ev_a120<-subset(ev_a120, ev_a120$effsize_0vs120_lower_CI < 0
                & ev_a120$effsize_0vs120_upper_CI < 0 |
                  ev_a120$effsize_0vs120_lower_CI > 0
                & ev_a120$effsize_0vs120_upper_CI > 0)

ev_a320<-rna_ev_t_test_and_effect_size_all_Sum
ev_a320<-subset(ev_a320, ev_a320$effsize_0vs320_lower_CI < 0
                & ev_a320$effsize_0vs320_upper_CI < 0 |
                  ev_a320$effsize_0vs320_lower_CI > 0
                & ev_a320$effsize_0vs320_upper_CI > 0)

View(rna_ev_t_test_and_effect_size_all_Sum)
View(ev_a120)
View(ev_a320)

#just in case I decide to use it in future, merge rnas present in both df's
ev_a120_320<-merge(ev_a120,ev_a320, by="Gene_ID", no.dups = TRUE)
View(ev_a120_320)

#you get duplicate columns, so solve that by choosing Proteins in a120 that are also present in a320:
ev_a120_320_1<-subset(ev_a120,ev_a120$Gene_ID %in%ev_a320$Gene_ID)
View(ev_a120_320_1)

getwd()

write.csv(ev_a120,file="rna_ev_t_test_and_effect_size_nonzero_CI_120.csv")
write.csv(ev_a320,file="rna_ev_t_test_and_effect_size_nonzero_CI_320.csv")
write.csv(ev_a120_320_1,file="rna_ev_t_test_and_effect_size_nonzero_CI_both_120_320.csv")


#cell now#
cell_a120<-rna_Cell_t_test_and_effect_size_all_Sum
cell_a120<-subset(cell_a120, cell_a120$effsize_0vs120_lower_CI < 0
                  & cell_a120$effsize_0vs120_upper_CI < 0 |
                    cell_a120$effsize_0vs120_lower_CI > 0
                  & cell_a120$effsize_0vs120_upper_CI > 0)

cell_a320<-rna_Cell_t_test_and_effect_size_all_Sum
cell_a320<-subset(cell_a320, cell_a320$effsize_0vs320_lower_CI < 0
                  & cell_a320$effsize_0vs320_upper_CI < 0 |
                    cell_a320$effsize_0vs320_lower_CI > 0
                  & cell_a320$effsize_0vs320_upper_CI > 0)

View(rna_Cell_t_test_and_effect_size_all_Sum)
View(cell_a120)
View(cell_a320)

#merge rnas present in both df's
cell_a120_320<-merge(cell_a120,cell_a320, by="Gene_ID", no.dups = TRUE)
View(cell_a120_320)

#you get duplicate columns, so solve that by choosing rnas in a120 that are also present in a320:
cell_a120_320_1<-subset(cell_a120,cell_a120$Gene_ID %in%cell_a320$Gene_ID)
View(cell_a120_320_1)

getwd()

write.csv(cell_a120,file="rna_Cell_t_test_and_effect_size_nonzero_CI_120.csv")
write.csv(cell_a320,file="rna_Cell_t_test_and_effect_size_nonzero_CI_320.csv")
write.csv(cell_a120_320_1,file="rna_Cell_t_test_and_effect_size_nonzero_CI_both_120_320.csv")

#for ev
length(which(ev_a120$pvalue_0vs120<0.05))
#671 out of 819 nonzero-containing 95% CI
length(which(ev_a120$pvalue_0vs120<0.05 & ev_a120$effsize_0vs120>0.4))
#305 has nonzero-containing 95% CI, pvalue<0.05, and effsize>0.04
length(which(ev_a120$pvalue_0vs120<0.05 & ev_a120$effsize_0vs120<(-0.4)))
#286 has nonzero-containing 95% CI, pvalue<0.05, and effsize<-0.04

length(which(ev_a320$pvalue_0vs320<0.05))
#1433 out of 1794 nonzeron-containing 95% CI
length(which(ev_a320$pvalue_0vs320<0.05 & ev_a320$effsize_0vs320>0.4))
#494 has nonzero-containing 95% CI, pvalue<0.05, and effsize>0.04
length(which(ev_a320$pvalue_0vs320<0.05 & ev_a320$effsize_0vs320<(-0.4)))
#556 has nonzero-containing 95% CI, pvalue<0.05, and effsize<-0.04

#for cell
length(which(cell_a120$pvalue_0vs120<0.05))
#852 out of 1119 nonzero-containing 95% CI
length(which(cell_a120$pvalue_0vs120<0.05 & cell_a120$effsize_0vs120>0.4))
#284 has nonzero-containing 95% CI, pvalue<0.05, and effsize>0.04
length(which(cell_a120$pvalue_0vs120<0.05 & cell_a120$effsize_0vs120<(-0.4)))
#306 has nonzero-containing 95% CI, pvalue<0.05, and effsize<-0.04

length(which(cell_a320$pvalue_0vs320<0.05))
#1646 out of 1986 nonzeron-containing 95% CI
length(which(cell_a320$pvalue_0vs320<0.05 & cell_a320$effsize_0vs320>0.4))
#144 has nonzero-containing 95% CI, pvalue<0.05, and effsize>0.04
length(which(cell_a320$pvalue_0vs320<0.05 & cell_a320$effsize_0vs320<(-0.4)))
#1171 has nonzero-containing 95% CI, pvalue<0.05, and effsize<-0.04

####
##let's see how many ethanol-sensitive RNAs are sex-biased DEGs in EVs and in cells.####
#remember that i wrote ev_a120 as rna_ev_t_test_and_effect_size_nonzero_CI_120.csv, etc,
#so keep that in mind that if i import these files back into R, it won't be ev_a120, 
#but rather, rna_ev_t_test_and_effect_size_nonzero_CI_120. 
#write.csv(ev_a120,file="rna_ev_t_test_and_effect_size_nonzero_CI_120.csv")
#write.csv(ev_a320,file="rna_ev_t_test_and_effect_size_nonzero_CI_320.csv")
#write.csv(ev_a120_320_1,file="rna_ev_t_test_and_effect_size_nonzero_CI_both_120_320.csv")

View(rna_ev_t_test_and_effect_size_nonzero_CI_120)
View(rna_ev_t_test_and_effect_size_nonzero_CI_320)
View(rna_Cell_t_test_and_effect_size_nonzero_CI_120)
View(rna_Cell_t_test_and_effect_size_nonzero_CI_320)

#################################################################################################
####lets do some ggplot of effsize vs expression level####
#find the median of normalized expression level
formatC(median(Protein_Sum_ev22),format = "e", digits = 2)
# "9.31e-04"
formatC(median(Protein_Sum_cell22),format = "e", digits = 2)
# new value  "1.47e-03"

#==============
# LOAD PACKAGES
#==============
install.packages("tidyverse")

library(ggplot2)
library(plyr)
library(tidyverse)

install.packages("plotly")
library(plotly)


library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions
install.packages("RColorBrewer")
library("RColorBrewer")

#colors http://www.sthda.com/english/wiki/colors-in-r


#now, previously when i used ggplot for volcano plots,
#some red dots are overlaying the white dots, while its the other way for other red dots.
#this is because the values are being plotted in order of the data. 
#if you want certain ones to be overlayed on top of others, you have to rearrange the data, 
#so r plots it in order, from top to bottom
#so, i created a d.f., rearranging by decreasing pvalue and plotted
#i also set x-axis limit to center the graph to 0 x-axis. 
rna_high_to_low_pvalue_120<-arrange(ev_a120,desc(pvalue_0vs120))
View(rna_high_to_low_pvalue_120)
#median of all rna_ev_Sum
formatC(median(rna_ev_Sum),format = "e", digits = 2)
#49.00342
#"4.90e+01"
View(rna_ev_Sum)
#median of just ev_a120, which is nonzero containing 95% CI
formatC(median(ev_a120$Total_Normalized_rna_Expression_Level_EV),format = "e", digits = 2)
#[1] "4.94e+02"

ev_a120<-rna_ev_t_test_and_effect_size_all_Sum
ev_a120<-subset(ev_a120, ev_a120$effsize_0vs120_lower_CI < 0
                & ev_a120$effsize_0vs120_upper_CI < 0 |
                  ev_a120$effsize_0vs120_lower_CI > 0
                & ev_a120$effsize_0vs120_upper_CI > 0)


require("ggrepel")
set.seed(09192022)     
options(ggrepel.max.overlaps = Inf)
rna_high_to_low_pvalue_120 %>% 
  mutate(P_not_significant = ifelse(pvalue_0vs120 > 0.05 & pvalue_0vs120 > 0.05,T, F)) %>%
  mutate(Significant = ifelse(pvalue_0vs120 < 0.05 & abs(effsize_0vs120) > 0.4,T, F)) %>% 
  ggplot(rna_high_to_low_pvalue_120 %>%
           aes(x = effsize_0vs120, y = `Total_Normalized_rna_Expression_Level_EV`,colour=Significant)) +
  geom_vline(xintercept=c(-0.8,-0.4,-0.2,0.2,0.4,0.8), linetype="dotted", color="violet", size=1) +
  geom_hline(yintercept=4.94e+02, linetype="dotted", color="blue", size=1.5) +
  geom_point(aes(color = Significant, fill = Significant), shape=21, size=2) +
  scale_fill_manual(values=c("white", "red")) + 
  scale_colour_manual(values = c('light grey', 'red')) +
  scale_x_continuous(name="Effect Size 0vs120",limits=c(-2.0,2.0),breaks=c(-2.0,-1.5,-1,-0.8,-0.4,0, 0.4,0.8, 1,1.5,2.0)) +
  scale_y_continuous(trans='log10') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits = c(10^0,10^5))+
  theme_bw()



getwd()
#i can control the y-axis if i want
#scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#labels = trans_format("log10", math_format(10^.x)), limits = c(10^-5,10^-0))

#to only have text for significant proteins
require("ggrepel")
set.seed(09192022)     
options(ggrepel.max.overlaps = Inf)
rna_high_to_low_pvalue_120 %>% 
  mutate(P_not_significant = ifelse(pvalue_0vs120 > 0.05 & pvalue_0vs120 > 0.05,T, F)) %>%
  mutate(Significant = ifelse(pvalue_0vs120 < 0.05 & abs(effsize_0vs120) > 0.4,T, F)) %>% 
  ggplot(rna_high_to_low_pvalue_120 %>%
           aes(x = effsize_0vs120, y = `Total_Normalized_rna_Expression_Level_EV`,colour=Significant)) +
  geom_vline(xintercept=c(-0.8,-0.4,-0.2,0.2,0.4,0.8), linetype="dotted", color="violet", size=1) +
  geom_hline(yintercept=4.94e+02, linetype="dotted", color="blue", size=1.5) +
  geom_point(aes(color = Significant, fill = Significant), shape=21, size=2) +
  geom_text_repel(aes(label=ifelse(pvalue_0vs120<0.05&abs(effsize_0vs120)>1.2,Gene_ID,"")))+
  scale_fill_manual(values=c("white", "red")) + 
  scale_colour_manual(values = c('light grey', 'red')) +
  scale_x_continuous(name="Effect Size 0vs120",limits=c(-2.0,2.0),breaks=c(-2.0,-1.5,-1,-0.8,-0.4,0, 0.4,0.8, 1,1.5,2.0)) +
  scale_y_continuous(trans='log10') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits = c(10^0,10^5)) +
  theme_bw()

#ev320
rna_high_to_low_pvalue_320<-arrange(ev_a320,desc(pvalue_0vs320))

#median of just ev_a320
formatC(median(ev_a320$Total_Normalized_rna_Expression_Level_EV),format = "e", digits = 2)
#[1] "2.61e+03"


require("ggrepel")
set.seed(09192022)     
options(ggrepel.max.overlaps = Inf)
rna_high_to_low_pvalue_320 %>% 
  mutate(P_not_significant = ifelse(pvalue_0vs320 > 0.05 & pvalue_0vs320 > 0.05,T, F)) %>%
  mutate(Significant = ifelse(pvalue_0vs320 < 0.05 & abs(effsize_0vs320) > 0.4,T, F)) %>% 
  ggplot(rna_high_to_low_pvalue_320 %>%
           aes(x = effsize_0vs320, y = `Total_Normalized_rna_Expression_Level_EV`,colour=Significant)) +
  geom_vline(xintercept=c(-0.8,-0.4,-0.2,0.2,0.4,0.8), linetype="dotted", color="violet", size=1) +
  geom_hline(yintercept=2.61e+03, linetype="dotted", color="blue", size=1.5) +
  geom_point(aes(color = Significant, fill = Significant), shape=21, size=2) +
  scale_fill_manual(values=c("white", "red")) + 
  scale_colour_manual(values = c('light grey', 'red')) +
  scale_x_continuous(name="Effect Size 0vs320",limits=c(-2.0,2.0),breaks=c(-2.0,-1.5,-1,-0.8,-0.4,0, 0.4,0.8, 1,1.5,2.0)) +
  scale_y_continuous(trans='log10') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits = c(10^0,10^5)) +
  theme_bw()



#### let's do one for cell
rna_cell_high_to_low_pvalue_120<-arrange(cell_a120,desc(pvalue_0vs120))

#median of all rna_ev_Sum
formatC(median(rna_Cell_Sum),format = "e", digits = 2)
#"2.13e+01"

#median of all 18 cell samples' rnas that are present just in cell_a120, meaning nonzero 95% CI
formatC(median(cell_a120$Total_Normalized_rna_Expression_Level_Cell),format = "e", digits = 2)
#[1] "1.97e+02"


require("ggrepel")
set.seed(09192022)     
options(ggrepel.max.overlaps = Inf)
rna_cell_high_to_low_pvalue_120 %>% 
  mutate(P_not_significant = ifelse(pvalue_0vs120 > 0.05 & pvalue_0vs120 > 0.05,T, F)) %>%
  mutate(Significant = ifelse(pvalue_0vs120 < 0.05 & abs(effsize_0vs120) > 0.4,T, F)) %>% 
  ggplot(rna_cell_high_to_low_pvalue_120 %>%
           aes(x = effsize_0vs120, y = `Total_Normalized_rna_Expression_Level_Cell`,colour=Significant)) +
  geom_vline(xintercept=c(-0.8,-0.4,-0.2,0.2,0.4,0.8), linetype="dotted", color="violet", size=1) +
  geom_hline(yintercept=1.97e+02, linetype="dotted", color="blue", size=1.5) +
  geom_point(aes(color = Significant, fill = Significant), shape=21, size=2) +
  scale_fill_manual(values=c("white", "red")) + 
  scale_colour_manual(values = c('light grey', 'red')) +
  scale_x_continuous(name="Effect Size 0vs120",limits=c(-2.0,2.0),breaks=c(-2.0,-1.5,-1,-0.8,-0.4,0, 0.4,0.8, 1,1.5,2.0)) +
  scale_y_continuous(trans='log10') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits = c(10^0,10^5)) +
  theme_bw()

#cell320
rna_cell_high_to_low_pvalue_320<-arrange(cell_a320,desc(pvalue_0vs320))

#median of all rna_ev_Sum
formatC(median(rna_Cell_Sum),format = "e", digits = 2)
#"2.13e+01"

#median of just ev_a120
formatC(median(cell_a320$Total_Normalized_rna_Expression_Level_Cell),format = "e", digits = 2)
#[1] "2.47e+03"

require("ggrepel")
set.seed(09192022)     
options(ggrepel.max.overlaps = Inf)
rna_cell_high_to_low_pvalue_320 %>% 
  mutate(P_not_significant = ifelse(pvalue_0vs320 > 0.05 & pvalue_0vs320 > 0.05,T, F)) %>%
  mutate(Significant = ifelse(pvalue_0vs320 < 0.05 & abs(effsize_0vs320) > 0.4,T, F)) %>% 
  ggplot(rna_cell_high_to_low_pvalue_320 %>%
           aes(x = effsize_0vs320, y = `Total_Normalized_rna_Expression_Level_Cell`,colour=Significant)) +
  geom_vline(xintercept=c(-0.8,-0.4,-0.2,0.2,0.4,0.8), linetype="dotted", color="violet", size=1) +
  geom_hline(yintercept=2.47e+03, linetype="dotted", color="blue", size=1.5) +
  geom_point(aes(color = Significant, fill = Significant), shape=21, size=2) +
  scale_fill_manual(values=c("white", "red")) + 
  scale_colour_manual(values = c('light grey', 'red')) +
  scale_x_continuous(name="Effect Size 0vs320",limits=c(-2.0,2.0),breaks=c(-2.0,-1.5,-1,-0.8,-0.4,0, 0.4,0.8, 1,1.5,2.0)) +
  scale_y_continuous(trans='log10') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)), limits = c(10^0,10^5)) +
  theme_bw()





###rna enrichment, looking at ratio of ev to cell and etoh effects on that enrichment####
###########
#let's do enrichment now.
#do a ratio of each rna's (ev/(ev+cell))
#let's do paired t test and effect size for enrichment data. 
#here, positive effect size means that ev enrichment has increased in treatment group compared to ctrl.
#it does not mean that the rna is enriched in Ev compared to cell.
#that you have to look at enrichment ratio, whether it is above 0.5
#but it does show whether the enrichment change in EV is an increase or a decrease of that rna in EV by treatment group.


#let's do some enrichment analysis
#the enrichment equation will be Enrichment = ev/(ev+cell)
#this way, i wont have any NaN or Inf
#if the ratio is 0.5>x>1, then it is enriched in ev's
#if the ratio is 0<x<0.5, then it is enriched in cells.

#lets do enrichment. Enrichment = ev/(ev+cell)
View(rna_expression_ev_cell)

rna_enrichment_ev_cell<-data.frame(matrix(NA, nrow = 18, ncol = ncol(rna_expression_ev_cell)))
for(i in 1:ncol(rna_expression_ev_cell)){rna_enrichment_ev_cell[,i]<-
  (rna_expression_ev_cell[1:18,i]/(rna_expression_ev_cell[1:18,i]+rna_expression_ev_cell[19:36,i]))}

colnames(rna_enrichment_ev_cell)<-colnames(rna_expression_ev_cell)[1:ncol(rna_expression_ev_cell)]
View(rna_enrichment_ev_cell)

rna_enrichment_ev_cell[,c(40206:40208)]<-enrichment_ntall_ev22_cell22[,c(1,2,3)]
rna_enrichment_ev_cell<-rna_enrichment_ev_cell[,c(40206:40208,1:40205)]
colnames(rna_enrichment_ev_cell)[2]<-"Set"
write.csv(rna_enrichment_ev_cell, file="rna_enrichment_ev_cell.csv")

##################################################################################

###as before, let's do p value and effsize of nonzero-containing 95% CI.
#because the sample size is small and variance large, most of the rnas have overlapping CI's when doing ctrl vs 120 or 320.
#for small samples, people do CI on the hedge's g for the effect size and see if it is nonzero standing 
#(anything that passes zero is no good; so ci's lower and upper limit need to have the same positive or negative values). 
#So I have two tests, 
#I am doing pvalue and effsize that includes nonzero standing CI from the Hedges's g. 
#So even though I am using two tests, I am strengthening the effsize test. 
#This is because my sample size is very small, 
#so the variance between samples is too wide and might not be the best test for my sample size. 
#also, we cannot do t.test on genes with NaN values (NaN was created when dividing something by 0, like 0/(0+0) for enrichment)
#thus, lets remove all genes (columns) that contain NaN for any sample (rows) first before doing t.test


rna_enrichment_ev_cell_no_nan<-
  rna_enrichment_ev_cell[,!colSums(is.na(rna_enrichment_ev_cell[1:ncol(rna_enrichment_ev_cell)]))]
View(rna_enrichment_ev_cell_no_nan)
#20922 genes left



library(effsize)

#ctrl vs 120

#practice
View(rna_enrichment_ev_cell)
ef_size<-cohen.d(rna_enrichment_ev_cell[c(2,5,8,11,14,17),7],rna_enrichment_ev_cell[c(1,4,7,10,13,16),7],pooled=TRUE,paired=TRUE,
                 na.rm=TRUE, hedges.correction=TRUE)
ef_size
View(ef_size)
#for effsize:
ef_size$estimate
#for effsize's lower conf interval:
ef_size$conf.int[[1]]
#for effsize's upper conf interval:
ef_size$conf.int[[2]]
#cool.
#let's do it for real:

# Convert all variable types to numeric
rna_enrichment_ev_cell_no_nan[,4:ncol(rna_enrichment_ev_cell_no_nan)] <- 
  as.data.frame(apply(rna_enrichment_ev_cell_no_nan[,4:ncol(rna_enrichment_ev_cell_no_nan)], 2, as.numeric)) 

tp_table_rna_enrichment_ctrl_vs_120<-data.frame(Gene_ID = rep(NA,ncol(rna_enrichment_ev_cell_no_nan)), 
                                                pvalue = rep(NA,ncol(rna_enrichment_ev_cell_no_nan)),
                                                effsize=rep(NA,ncol(rna_enrichment_ev_cell_no_nan)),
                                                effsize_lower_CI=rep(NA,ncol(rna_enrichment_ev_cell_no_nan)),
                                                effsize_upper_CI=rep(NA,ncol(rna_enrichment_ev_cell_no_nan)),
                                                mean_control=rep(NA,ncol(rna_enrichment_ev_cell_no_nan)),
                                                mean_ethanol=rep(NA,ncol(rna_enrichment_ev_cell_no_nan)))


#rna starts from column 4 in this case, hence i in 4:
for (i in 4:ncol(rna_enrichment_ev_cell_no_nan))
{evttest<-t.test(rna_enrichment_ev_cell_no_nan[c(1,4,7,10,13,16),i],
                 rna_enrichment_ev_cell_no_nan[c(2,5,8,11,14,17),i],
                 paired=TRUE)
tp_table_rna_enrichment_ctrl_vs_120[i,2]<-evttest$p.value
tp_table_rna_enrichment_ctrl_vs_120[i,1]<-colnames(rna_enrichment_ev_cell_no_nan)[i]
ef_size<-cohen.d(rna_enrichment_ev_cell_no_nan[c(2,5,8,11,14,17),i],
                 rna_enrichment_ev_cell_no_nan[c(1,4,7,10,13,16),i],
                 pooled=TRUE,paired=TRUE,na.rm=TRUE, hedges.correction=TRUE)
tp_table_rna_enrichment_ctrl_vs_120[i,3]<-ef_size$estimate
tp_table_rna_enrichment_ctrl_vs_120[i,4]<-ef_size$conf.int[[1]]
tp_table_rna_enrichment_ctrl_vs_120[i,5]<-ef_size$conf.int[[2]]
tp_table_rna_enrichment_ctrl_vs_120[i,7]<-mean(rna_enrichment_ev_cell_no_nan[c(2,5,8,11,14,17),i])
tp_table_rna_enrichment_ctrl_vs_120[i,6]<-mean(rna_enrichment_ev_cell_no_nan[c(1,4,7,10,13,16),i])}
View(tp_table_rna_enrichment_ctrl_vs_120)


#now 0vs320
tp_table_rna_enrichment_ctrl_vs_320<-data.frame(Gene_ID = rep(NA,ncol(rna_enrichment_ev_cell_no_nan)), 
                                                pvalue = rep(NA,ncol(rna_enrichment_ev_cell_no_nan)),
                                                effsize=rep(NA,ncol(rna_enrichment_ev_cell_no_nan)),
                                                effsize_lower_CI=rep(NA,ncol(rna_enrichment_ev_cell_no_nan)),
                                                effsize_upper_CI=rep(NA,ncol(rna_enrichment_ev_cell_no_nan)),
                                                mean_control=rep(NA,ncol(rna_enrichment_ev_cell_no_nan)),
                                                mean_ethanol=rep(NA,ncol(rna_enrichment_ev_cell_no_nan)))

#rna starts from column 1 in this case, hence i in 1:
for (i in 4:ncol(rna_enrichment_ev_cell_no_nan))
{evttest<-t.test(rna_enrichment_ev_cell_no_nan[c(1,4,7,10,13,16),i],
                 rna_enrichment_ev_cell_no_nan[c(3,6,9,12,15,18),i],
                 paired=TRUE)
tp_table_rna_enrichment_ctrl_vs_320[i,2]<-evttest$p.value
tp_table_rna_enrichment_ctrl_vs_320[i,1]<-colnames(rna_enrichment_ev_cell_no_nan)[i]
ef_size<-cohen.d(rna_enrichment_ev_cell_no_nan[c(3,6,9,12,15,18),i],
                 rna_enrichment_ev_cell_no_nan[c(1,4,7,10,13,16),i],
                 pooled=TRUE,paired=TRUE,na.rm=TRUE, hedges.correction=TRUE)
tp_table_rna_enrichment_ctrl_vs_320[i,3]<-ef_size$estimate
tp_table_rna_enrichment_ctrl_vs_320[i,4]<-ef_size$conf.int[[1]]
tp_table_rna_enrichment_ctrl_vs_320[i,5]<-ef_size$conf.int[[2]]
tp_table_rna_enrichment_ctrl_vs_320[i,7]<-mean(rna_enrichment_ev_cell_no_nan[c(3,6,9,12,15,18),i])
tp_table_rna_enrichment_ctrl_vs_320[i,6]<-mean(rna_enrichment_ev_cell_no_nan[c(1,4,7,10,13,16),i])}
dim(tp_table_rna_enrichment_ctrl_vs_320)
View(tp_table_rna_enrichment_ctrl_vs_320)
summary(tp_table_rna_enrichment_ctrl_vs_320$mean_control)
summary(tp_table_rna_enrichment_ctrl_vs_320$mean_ethanol)

#let's see how i can merge the data, 
#so I can have 0vs320 pvalue and effect size next to 0vs120 proteins with significant pvalue, and vice versa.
#better yet, lets merge the two d.f. together.
#to do that, let's rename some columns to specify which treatment it comes from. 
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))

library(data.table)

setnames(tp_table_rna_enrichment_ctrl_vs_120, old=c("pvalue","effsize","effsize_lower_CI","effsize_upper_CI","mean_control","mean_ethanol"), 
         new=c("pvalue_0vs120", "effsize_0vs120","effsize_0vs120_lower_CI","effsize_0vs120_upper_CI","mean_control","mean_ethanol_120"))
setnames(tp_table_rna_enrichment_ctrl_vs_320, old=c("pvalue","effsize","effsize_lower_CI","effsize_upper_CI","mean_control","mean_ethanol"), 
         new=c("pvalue_0vs320", "effsize_0vs320","effsize_0vs320_lower_CI","effsize_0vs320_upper_CI","mean_control","mean_ethanol_320"))
#get rid of empty 3 rows
tp_table_rna_enrichment_ctrl_vs_120<-tp_table_rna_enrichment_ctrl_vs_120[-c(1:3),]
tp_table_rna_enrichment_ctrl_vs_320<-tp_table_rna_enrichment_ctrl_vs_320[-c(1:3),]

#for merging into a new d.f., let's only include columns I am interested in merging into. 
rna_enrichment_t_test_and_effect_size_all<-merge(tp_table_rna_enrichment_ctrl_vs_120,tp_table_rna_enrichment_ctrl_vs_320,by="Gene_ID")
View(rna_enrichment_t_test_and_effect_size_all)

#there is mean_ctrl.x and mean_ctrl.y with same values so delete one
rna_enrichment_t_test_and_effect_size_all<-rna_enrichment_t_test_and_effect_size_all[,-c(12)]

# Rename a column in R
colnames(rna_enrichment_t_test_and_effect_size_all)[colnames(rna_enrichment_t_test_and_effect_size_all)=="mean_control.x"] <- "mean_ctrl"


#remove duplicates
#it worked. but there are duplicate rows, so let's remove that with the function duplicated() and unique() to extract only unique rows
library(tidyverse)
#If you want to remove duplicated elements, use !duplicated(), where ! is a logical negation, telling we don't want duplicate rows.
#Remove duplicates based on Accession Number column, since R doesnt recognize space, rename Accession Number to Accession_Number
#colnames(rna_ev_t_test_and_effect_size_all_Sum)[3]<-"Accession_Number"
rna_enrichment_t_test_and_effect_size_all<-rna_enrichment_t_test_and_effect_size_all %>% distinct(Gene_ID, .keep_all = TRUE)
View(rna_enrichment_t_test_and_effect_size_all)

write.csv(rna_enrichment_t_test_and_effect_size_all,file="rna_enrichment_t_test_and_effect_size_all.csv")


#now, let's create a d.f. with only nonzero CI rnas.
#we can use subset() with multiple conditions, with & for "and", and | for "or".

a120<-rna_enrichment_t_test_and_effect_size_all

a120<-subset(a120, a120$effsize_0vs120_lower_CI < 0
             & a120$effsize_0vs120_upper_CI < 0 |
               a120$effsize_0vs120_lower_CI > 0
             & a120$effsize_0vs120_upper_CI > 0)

a320<-rna_enrichment_t_test_and_effect_size_all
a320<-subset(a320, a320$effsize_0vs320_lower_CI < 0
             & a320$effsize_0vs320_upper_CI < 0 |
               a320$effsize_0vs320_lower_CI > 0
             & a320$effsize_0vs320_upper_CI > 0)

View(a120)
View(a320)

length(which(a120$pvalue_0vs120<0.05))
#632 out of 763 nonzero-containing 95% CI
length(which(a120$pvalue_0vs120<0.05 & a120$effsize_0vs120>0.4))
#278 has nonzero-containing 95% CI, pvalue<0.05, and effsize>0.04
length(which(a120$pvalue_0vs120<0.05 & a120$effsize_0vs120<(-0.4)))
#241 has nonzero-containing 95% CI, pvalue<0.05, and effsize<-0.04

length(which(a320$pvalue_0vs320<0.05))
#926 out of 1129 nonzeron-containing 95% CI
length(which(a320$pvalue_0vs320<0.05 & a320$effsize_0vs320>0.4))
#538 has nonzero-containing 95% CI, pvalue<0.05, and effsize>0.04
length(which(a320$pvalue_0vs320<0.05 & a320$effsize_0vs320<(-0.4)))
#196 has nonzero-containing 95% CI, pvalue<0.05, and effsize<-0.04

#merge proteins present in both df's
a120_320<-merge(a120,a320, by="Gene_ID", no.dups = TRUE)
View(a120_320)

#you get duplicate columns, so solve that by choosing Proteins in a120 that are also present in a320:
a120_320_1<-subset(a120,a120$Gene_ID %in%a320$Gene_ID)
View(a120_320_1)

getwd()

#let's see how many ethanol-sensitive RNAs are also sex-dependent DEGs in EVs or in Cells
ev_120
View(rna_ev_t_test_and_effect_size_nonzero_CI_120)
View(rna_ev_t_test_and_effect_size_nonzero_CI_320)
View(rna_Cell_t_test_and_effect_size_nonzero_CI_120)
View(rna_Cell_t_test_and_effect_size_nonzero_CI_320)

ev120_p0.05_nonzero_effsize0.4<-rna_ev_t_test_and_effect_size_nonzero_CI_120
ev320_p0.05_nonzero_effsize0.4<-rna_ev_t_test_and_effect_size_nonzero_CI_320
Cell120_p0.05_nonzero_effsize0.4<-rna_Cell_t_test_and_effect_size_nonzero_CI_120
Cell320_p0.05_nonzero_effsize0.4<-rna_Cell_t_test_and_effect_size_nonzero_CI_320

ev120_p0.05_nonzero_effsize0.4<-
  subset(ev120_p0.05_nonzero_effsize0.4, ev120_p0.05_nonzero_effsize0.4$pvalue_0vs120<0.05
         & abs(ev120_p0.05_nonzero_effsize0.4$effsize_0vs120) > 0.4)
ev320_p0.05_nonzero_effsize0.4<-
  subset(ev320_p0.05_nonzero_effsize0.4, ev320_p0.05_nonzero_effsize0.4$pvalue_0vs320<0.05
         & abs(ev320_p0.05_nonzero_effsize0.4$effsize_0vs320) > 0.4)

Cell120_p0.05_nonzero_effsize0.4<-
  subset(Cell120_p0.05_nonzero_effsize0.4, Cell120_p0.05_nonzero_effsize0.4$pvalue_0vs120<0.05
         & abs(Cell120_p0.05_nonzero_effsize0.4$effsize_0vs120) > 0.4)
Cell320_p0.05_nonzero_effsize0.4<-
  subset(Cell320_p0.05_nonzero_effsize0.4, Cell320_p0.05_nonzero_effsize0.4$pvalue_0vs320<0.05
         & abs(Cell320_p0.05_nonzero_effsize0.4$effsize_0vs320) > 0.4)

View(ev120_p0.05_nonzero_effsize0.4)
View(ev320_p0.05_nonzero_effsize0.4)
View(Cell120_p0.05_nonzero_effsize0.4)
View(Cell320_p0.05_nonzero_effsize0.4)

ev120_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV<-
  subset(ev120_p0.05_nonzero_effsize0.4,ev120_p0.05_nonzero_effsize0.4$Gene_ID %in%DESEQ2_results_EV_sig_genes$Gene_ID)
View(ev120_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV)
#From the 591 moderate ethanol-sensitive RNAs in EVs, 55 of them are significant DEGs by sex in EVs
ev120_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell<-
  subset(ev120_p0.05_nonzero_effsize0.4,ev120_p0.05_nonzero_effsize0.4$Gene_ID %in%DESEQ2_results_Cell_sig_genes$Gene_ID)
View(ev120_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell)
#From the 591 moderate ethanol-sensitive RNAs in EVs, 48 of them are significant DEGs by sex in cells

ev320_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV<-
  subset(ev320_p0.05_nonzero_effsize0.4,ev320_p0.05_nonzero_effsize0.4$Gene_ID %in%DESEQ2_results_EV_sig_genes$Gene_ID)
View(ev320_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV)
#From the 1050 heavy ethanol-sensitive RNAs in EVs, 118 of them are significant DEGs by sex in EVs
ev320_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell<-
  subset(ev320_p0.05_nonzero_effsize0.4,ev320_p0.05_nonzero_effsize0.4$Gene_ID %in%DESEQ2_results_Cell_sig_genes$Gene_ID)
View(ev320_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell)
#From the 1050 heavy ethanol-sensitive RNAs in EVs, 194 of them are significant DEGs by sex in cells

Cell120_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell<-
  subset(Cell120_p0.05_nonzero_effsize0.4,Cell120_p0.05_nonzero_effsize0.4$Gene_ID %in%DESEQ2_results_Cell_sig_genes$Gene_ID)
View(Cell120_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell)
#From the 590 moderate ethanol-sensitive RNAs in Cells, 48 of them are significant DEGs by sex in Cells

Cell320_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell<-
  subset(Cell320_p0.05_nonzero_effsize0.4,Cell320_p0.05_nonzero_effsize0.4$Gene_ID %in%DESEQ2_results_Cell_sig_genes$Gene_ID)
View(Cell320_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell)
#From the 1315 heavy ethanol-sensitive RNAs in Cells, 169 of them are significant DEGs by sex in Cells





View(rna_enrichment_t_test_and_effect_size_nonzero_CI_120)
View(rna_enrichment_t_test_and_effect_size_nonzero_CI_320)


length(which(rna_enrichment_t_test_and_effect_size_nonzero_CI_120$pvalue_0vs120<0.05))
#632 out of 763 nonzero-containing 95% CI
length(which(rna_enrichment_t_test_and_effect_size_nonzero_CI_120$pvalue_0vs120<0.05 & rna_enrichment_t_test_and_effect_size_nonzero_CI_120$effsize_0vs120>0.4))
#278 has nonzero-containing 95% CI, pvalue<0.05, and effsize>0.4
length(which(rna_enrichment_t_test_and_effect_size_nonzero_CI_120$pvalue_0vs120<0.05 & rna_enrichment_t_test_and_effect_size_nonzero_CI_120$effsize_0vs120<(-0.4)))
#241 has nonzero-containing 95% CI, pvalue<0.05, and effsize<-0.4

a120_p0.05_nonzero_effsize0.4<-rna_enrichment_t_test_and_effect_size_nonzero_CI_120
a120_p0.05_nonzero_effsize0.4<-
  subset(a120_p0.05_nonzero_effsize0.4, a120_p0.05_nonzero_effsize0.4$pvalue_0vs120<0.05
         & abs(a120_p0.05_nonzero_effsize0.4$effsize_0vs120) > 0.4)

a320_p0.05_nonzero_effsize0.4<-rna_enrichment_t_test_and_effect_size_nonzero_CI_320
a320_p0.05_nonzero_effsize0.4<-
  subset(a320_p0.05_nonzero_effsize0.4, a320_p0.05_nonzero_effsize0.4$pvalue_0vs320<0.05
         & abs(a320_p0.05_nonzero_effsize0.4$effsize_0vs320) > 0.4)

DESEQ2_results_EV_sig_genes
colnames(DESEQ2_results_EV_sig_genes)[1] <- "Gene_ID"
colnames(DESEQ2_results_Cell_sig_genes)[1] <- "Gene_ID"
a120_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV<-
  subset(a120_p0.05_nonzero_effsize0.4,a120_p0.05_nonzero_effsize0.4$Gene_ID %in%DESEQ2_results_EV_sig_genes$Gene_ID)
View(a120_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV)
#From the 519 moderate ethanol-sensitive RNAs for EV to cell ratio, 46 of them are significant DEGs by sex in EVs
a120_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell<-
  subset(a120_p0.05_nonzero_effsize0.4,a120_p0.05_nonzero_effsize0.4$Gene_ID %in%DESEQ2_results_Cell_sig_genes$Gene_ID)
View(a120_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell)
#From the 519 moderate ethanol-sensitive RNAs for EV to cell ratio, 41 of them are significant DEGs by sex in cells

a320_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV<-
  subset(a320_p0.05_nonzero_effsize0.4,a320_p0.05_nonzero_effsize0.4$Gene_ID %in%DESEQ2_results_EV_sig_genes$Gene_ID)
View(a320_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV)
#From the 734 heavy ethanol-sensitive RNAs for EV to cell ratio, 82 of them are significant DEGs by sex in EVs
a320_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell<-
  subset(a320_p0.05_nonzero_effsize0.4,a320_p0.05_nonzero_effsize0.4$Gene_ID %in%DESEQ2_results_Cell_sig_genes$Gene_ID)
View(a320_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell)
#From the 734 heavy ethanol-sensitive RNAs for EV to cell ratio, 95 of them are significant DEGs by sex in cells

#let's actually look at the ones that are enriched in EVs and depleted in cells,
#and see how many of them are sex-invariant and sex-variant response to alcohol
a120_p0.05_nonzero_effsize_greater_0.4<-
  subset(a120_p0.05_nonzero_effsize0.4, a120_p0.05_nonzero_effsize0.4$pvalue_0vs120<0.05
         & a120_p0.05_nonzero_effsize0.4$effsize_0vs120 > 0.4)

a320_p0.05_nonzero_effsize_greater_0.4<-
  subset(a320_p0.05_nonzero_effsize0.4, a320_p0.05_nonzero_effsize0.4$pvalue_0vs320<0.05
         & a320_p0.05_nonzero_effsize0.4$effsize_0vs320 > 0.4)

View(a120_p0.05_nonzero_effsize_greater_0.4)
View(a320_p0.05_nonzero_effsize_greater_0.4)


#From the 278 moderate ethanol-sensitive RNAs enriched in EVs and depleted in cells (EV:Cell), 
#only 19 of them are significant DEGs by sex in EVs

a320_p0.05_nonzero_effsize_greater_0.4_sig_DEGs_in_EV<-
  subset(a320_p0.05_nonzero_effsize_greater_0.4,a320_p0.05_nonzero_effsize_greater_0.4$Gene_ID %in%DESEQ2_results_EV_sig_genes$Gene_ID)
View(a320_p0.05_nonzero_effsize_greater_0.4_sig_DEGs_in_EV)
#From the 538 heavy ethanol-sensitive RNAs enriched in EVs and depleted in cells (EV:Cell), 
#only 57 of them are significant DEGs by sex in EVs



View(ev120_p0.05_nonzero_effsize0.4)
View(ev320_p0.05_nonzero_effsize0.4)
View(Cell120_p0.05_nonzero_effsize0.4)
View(Cell320_p0.05_nonzero_effsize0.4)

write.csv(ev120_p0.05_nonzero_effsize0.4,file="ev120_p0.05_nonzero_effsize0.4.csv")
write.csv(ev320_p0.05_nonzero_effsize0.4,file="ev320_p0.05_nonzero_effsize0.4.csv")
write.csv(Cell120_p0.05_nonzero_effsize0.4,file="Cell120_p0.05_nonzero_effsize0.4.csv")
write.csv(Cell320_p0.05_nonzero_effsize0.4,file="Cell320_p0.05_nonzero_effsize0.4.csv")
write.csv(a120_p0.05_nonzero_effsize0.4,file="EVtoCellenrichment120_p0.05_nonzero_effsize0.4.csv")
write.csv(a320_p0.05_nonzero_effsize0.4,file="EVtoCellenrichment320_p0.05_nonzero_effsize0.4.csv")

write.csv(ev120_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV,file="ev120_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV.csv")
write.csv(ev320_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV,file="ev320_p0.05_nonzero_effsize0.4_sig_DEGs_in_EV.csv")
write.csv(Cell120_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell,file="Cell120_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell.csv")
write.csv(Cell320_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell,file="Cell320_p0.05_nonzero_effsize0.4_sig_DEGs_in_Cell.csv")


write.csv(a120_p0.05_nonzero_effsize_greater_0.4_sig_DEGs_in_EV,file="EVtoCellenrichment120_p0.05_nonzero_effsize_greater_0.4_sex_variant.csv")
write.csv(a320_p0.05_nonzero_effsize_greater_0.4_sig_DEGs_in_EV,file="EVtoCellenrichment320_p0.05_nonzero_effsize_greater_0.4_sex_variant.csv")

write.csv(a120_p0.05_nonzero_effsize_greater_0.4,file="EVtoCellenrichment120_p0.05_nonzero_effsize_greater_0.4.csv")
write.csv(a320_p0.05_nonzero_effsize_greater_0.4,file="EVtoCellenrichment320_p0.05_nonzero_effsize_greater_0.4.csv")



write.csv(a120,file="rna_enrichment_t_test_and_effect_size_nonzero_CI_120.csv")
write.csv(a320,file="rna_enrichment_t_test_and_effect_size_nonzero_CI_320.csv")
write.csv(a120_320_1,file="rna_enrichment_t_test_and_effect_size_nonzero_CI_both_120_320.csv")

#choose sig pvalue in excel and do graphpad
#worked well, but now, we want to include those rnas that did not pass nonzero CI test.
#i will have 3 different colors: 1. red for psig. 2. blue for psig and nonzero CI. 3. gray or black for the rest. 

#create a df that is a subset of rnas that are not in a120, meaning no nonzero CI rnas.
a120_zero_CI<-rna_enrichment_t_test_and_effect_size_all
a120_zero_CI<-subset(a120_zero_CI,!a120_zero_CI$Gene_ID %in%a120$Gene_ID)

a320_zero_CI<-rna_enrichment_t_test_and_effect_size_all
a320_zero_CI<-subset(a320_zero_CI,!a320_zero_CI$Gene_ID %in%a320$Gene_ID)
View(a320_zero_CI)
write.csv(a120_zero_CI,file="rna_enrichment_t_test_and_effect_size_zero_CI_120.csv")
write.csv(a320_zero_CI,file="rna_enrichment_t_test_and_effect_size_zero_CI_320.csv")


getwd()
setwd("C:/Users/daehy/Dropbox/Transcriptomic/R")


######################################3

#let's do EnhancedVolcano
#y axis will be -log(pvalue) and x axis will be effsize
#first, row names from gene id
rownames(a120)<-a120$Gene_ID
#custom cut-offs: FC_log2 0.5-1.5; default adj_pval of 10e-6
#ev/cell 120vs0 first
set.seed(09192022)     
pdf("RNA_EV_Enrichment_120vs0_Volcano_plot_custom.pdf")
require("ggrepel")
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(a120,
                lab = rownames(a120),
                x = 'effsize_0vs120',
                y = 'pvalue_0vs120',
                xlim = c(-2, 2),
                ylim = c(0, 6),
                title = 'RNA Enrichment 120 vs 0 in EV to Cell',
                xlab = 'Effect Size 120 vs 0',
                pCutoff = 0.05,
                FCcutoff = 0.4,
                
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                legendLabels = c("NS","Effect Size","p-value","p-value and Effect Size"),
                colConnectors = 'grey30')
dev.off()

#or png
set.seed(09192022)     
png("RNA_EV_Enrichment_120vs0_Volcano_plot_custom.png", width = 2000, height = 2000, res = 300)

require("ggrepel")
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(a120,
                lab = rownames(a120),
                x = 'effsize_0vs120',
                y = 'pvalue_0vs120',
                xlim = c(-2, 2),
                ylim = c(0, 6),
                title = 'RNA Enrichment 120 vs 0 in EV to Cell',
                xlab = 'Effect Size 120 vs 0',
                pCutoff = 0.05,
                FCcutoff = 0.4,
                
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                legendLabels = c("NS","Effect Size","p-value","p-value and Effect Size"),
                colConnectors = 'grey30')
dev.off()

#ev/cell 320vs0 next
rownames(a320)<-a320$Gene_ID

#png

require("ggrepel")
set.seed(09192022)     
options(ggrepel.max.overlaps = Inf)
png("RNA_EV_Enrichment_320vs0_Volcano_plot_custom.png", width = 2000, height = 2000, res = 300)
EnhancedVolcano(a320,
                lab = rownames(a320),
                x = 'effsize_0vs320',
                y = 'pvalue_0vs320',
                xlim = c(-2, 2),
                ylim = c(0, 6),
                title = 'RNA Enrichment 320 vs 0 in EV to Cell',
                xlab = 'Effect Size 320 vs 0',
                pCutoff = 0.05,
                FCcutoff = 0.4,
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                legendLabels = c("NS","Effect Size","p-value","p-value and Effect Size"),
                colConnectors = 'grey30')
dev.off()

#or pdf
require("ggrepel")
set.seed(09192022)     
options(ggrepel.max.overlaps = Inf)
pdf("RNA_EV_Enrichment_320vs0_Volcano_plot_custom.pdf")
EnhancedVolcano(a320,
                lab = rownames(a320),
                x = 'effsize_0vs320',
                y = 'pvalue_0vs320',
                xlim = c(-2, 2),
                ylim = c(0, 6),
                title = 'RNA Enrichment 320 vs 0 in EV to Cell',
                xlab = 'Effect Size 320 vs 0',
                pCutoff = 0.05,
                FCcutoff = 0.4,
                
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                legendLabels = c("NS","Effect Size","p-value","p-value and Effect Size"),
                colConnectors = 'grey30')
dev.off()







####################################################################

####let's do pathway analysis. ####
####practice on default genelist if you want ####
#http://bioconductor.org/packages/release/bioc/html/ReactomePA.html
#https://guangchuangyu.github.io/software/ReactomePA/
#https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html#supported-organisms

#Reactome enrichment analysis
#ReactomePA is designed for reactome pathway based analysis (Yu and He 2016). 
#Reactome is an open-source, open access, manually curated and peer-reviewed pathway database.

#To install this package, start R (version "4.1") and enter:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")
library(ReactomePA)
#For older versions of R, please refer to the appropriate Bioconductor release.
#worked.

#Documentation
#To view documentation for the version of this package installed in your system, start R and enter:
browseVignettes("ReactomePA")


#We can translate from one type to other types.
#enrichment tool for interpreting omics data


#will need to convert to Enterz ID
#https://yulab-smu.github.io/clusterProfiler-book/chapter14.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
#install db for mouse IDs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DOSE")
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(DOSE)

library(ReactomePA)
#create vector of IDs
#a120 is rna_enrichment_t_test_and_effect_size_nonzero_CI_120
View(a120)
ID120 <- a120$Gene_ID
eg <- bitr(ID120, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg)
#  SYMBOL        ENTREZID
#1 0610006L08Rik     76253
#2 0610009B22Rik     66050
#3 0610009E02Rik 100125929
#4 0610009L18Rik     66838
#5 0610010F05Rik     71675
#6 0610010K14Rik    104457

#merge a120$effsize_0vs120 and eg
#these are rnas with nonzero CI for effectsize. 
#let's also make it where it will be these rnas but limited even more with psig value < 0.05
a120_psig<-a120
a120_psig<-subset(a120,!a120$pvalue_0vs120>0.05)

View(a120_psig)

View(a320_psig)
View(a320_Entrez_effsize)
View(a120_Entrez_effsize)
#rename column in eg so matching title
View(eg)
colnames(eg)[1] <- "Gene_ID"
#merge based on Gene_ID for effsize 0V120
a120_Entrez_effsize <- merge(eg, a120_psig[,c(1,3)], by = "Gene_ID")
View(a120_Entrez_effsize)
#get rid of accession number so i have entrez and effsize for 2 columns only
a120_Entrez_effsize <- a120_Entrez_effsize[,-1]
#make it numeric like geneList
#https://bioinformatics.stackexchange.com/questions/13036/reactomepa-dataset
myGeneList <- a120_Entrez_effsize[,2]
names(myGeneList) <- as.character(a120_Entrez_effsize[,1])
myGeneList <- sort(myGeneList, decreasing = TRUE)
head(myGeneList)
View(myGeneList)


#start processing for Reactome
#limiting proteins with effectsize that are greater than 0.5.
#use abs() if you dont care if the effectsize is negative or positive.
#class(geneList) #its numeric
#de <- names(myGeneList)[abs(myGeneList) > 0.5]
#head(de)
class(myGeneList) #its numeric
de <- names(myGeneList)[myGeneList > 0.4]
head(de)
View(as.data.frame(de))

x <- enrichPathway(gene=de,organism = "mouse", pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T)
head(as.data.frame(x))
#getting zero enriched pathways even when i increase the cutoff for effsize or pvalue. 
#this is a valid result.


###let's do one for a320##
#create vector of IDs
View(a320)

ID320 <- a320$Gene_ID
eg320 <- bitr(ID320, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg320)
View(eg320)
#these are rnas with nonzero CI for effectsize. 
#let's also make it where it will be these rnas but limited even more with psig value < 0.05
a320_psig<-subset(a320,!a320$pvalue_0vs320>0.05)
View(a320_psig)
View(a320_Entrez_effsize)

#merge a320$effsize_0vs320 (column 9 in this d.f.) and eg
#rename column in eg so matching title
colnames(eg320)[1] <- "Gene_ID"
#merge based on Accession_Number for effsize 0V320
a320_Entrez_effsize <- merge(eg320, a320_psig[,c(1,9)], by = "Gene_ID")
a320_Entrez_effsize <- a320_Entrez_effsize[,-1]
#make it numeric like geneList
#https://bioinformatics.stackexchange.com/questions/13036/reactomepa-dataset
myGeneList320 <- a320_Entrez_effsize[,2]
names(myGeneList320) <- as.character(a320_Entrez_effsize[,1])
myGeneList320 <- sort(myGeneList320, decreasing = TRUE)
head(myGeneList320)
View(myGeneList320)



#start processing for Reactome
class(myGeneList320) #its numeric
de320 <- names(myGeneList320)[(myGeneList320) > 0.4]
head(de320)
View(as.data.frame(de320))


x320 <- enrichPathway(gene=de320,organism = "mouse", pvalueCutoff=0.05, pAdjustMethod = "BH", readable=T)
head(as.data.frame(x320))
er320 <- as.data.frame(x320)
View(er320)
write.csv(er320, "rna_enrichment_a320_psig_effsize+0.4__nonzeroCI_Reactome_Results.csv")

png("rna_a320_bar_Reactome.png", width = 4000, height = 2000, res = 300)
barplot(x320, showCategory=10)
dev.off()

png("a320_dotplot_Reactome.png", width = 1700, height = 2000, res = 300)
dotplot(x320, showCategory=10)
dev.off()

png("a320_dotplot_Reactome.png", width = 2000, height = 2300, res = 300)
dotplot(x320, showCategory=13, font.size = 9)
dev.off()

RStudio.Version()

##https://guangchuangyu.github.io/2015/06/dotplot-for-enrichment-result/
##additional dot plot adjustments!!!!!!

x320_2 <- pairwise_termsim(x320) 
emapplot(x320_2, color="pvalue")

png("a320_enrichMap_Reactome.png", width = 2000, height = 2000, res = 300)
emapplot(x320_2)
dev.off()

png("a320_cnetPlot_Reactome.png", width = 2000, height = 2000, res = 300)
cnetplot(x320_2, categorySize="pvalue", foldChange=myGeneList320,
         cex_label_category = .8, 
         cex_label_gene = 1)
dev.off()

#GSEA
#Gene Set Enrichment Analysis (GSEA) is a computational method 
#that determines whether an a priori defined set of genes shows statistically
#significant, concordant differences between two biological states (e.g. phenotypes).
#y320 <- gsePathway(myGeneList320, nPerm=10000,
#                   organism = "mouse",
#                   pvalueCutoff=0.2,
#                   pAdjustMethod="BH", verbose=FALSE)

y320 <- gsePathway(myGeneList320, 
                   organism = "mouse",
                   pvalueCutoff=0.05,
                   pAdjustMethod="BH", verbose=FALSE)
#no term enriched under specific pvalueCutoff...

y320_2 <- pairwise_termsim(y320) 

png("rna_a320_GSEA_Reactome.png", width = 2000, height = 2000, res = 300)
emapplot(y320_2)
dev.off()


#view Specific pathway if you want.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("graphite")
library(graphite)


png("a320_RHO_GTPase_Effectors_Reactome.png", width = 2000, height = 2000, res = 300)
viewPathway("RHO GTPase Effectors", organism = "mouse", readable=TRUE, guide = "none", foldChange=myGeneList320)
dev.off()

png("a120_DNAReplicationPre-Initiation_Reactome.png", width = 2000, height = 2000, res = 300)
viewPathway("Axon guidance", organism = "mouse", readable=TRUE, guide = "none", foldChange=myGeneList320)
dev.off()

png("Hip_RecyclingL1_Reactome.png", width = 2000, height = 2000, res = 300)
viewPathway("Recycling pathway of L1", organism = "mouse", readable=TRUE, foldChange=myGeneList320)
dev.off()

########
#KEGG enrichment analysis
#visualize KEGG
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")
library("pathview")

#7.2 KEGG pathway over-representation analysis


de <- names(myGeneList)[(myGeneList) > 0.4]
kk <- enrichKEGG(gene         = de,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
head(kk)

View(as.data.frame(kk))
er_kk<- as.data.frame(kk)
write.csv(er_kk, "a120_psig_effsize+0.4_KEGG_Results.csv")

browseKEGG(kk, 'mmu00030')

#320
de320 <- names(myGeneList320)[(myGeneList320) > 0.4]
kk320 <- enrichKEGG(gene         = de320,
                    organism     = 'mmu',
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")
head(kk320)

View(as.data.frame(kk320))
er_kk320 <- as.data.frame(kk320)
write.csv(er_kk320, "a320_psig_effsize+0.4_KEGG_Results.csv")



#############################################
###WGCNA for transcriptomic paper####

save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")
View(datExpr)
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

#Weighted correlation network analysis (WGCNA) can be used for finding clusters (modules) of highly correlated genes, 
#for summarizing such clusters using the module eigengene or an intramodular hub gene, 
#for relating modules to one another and to external sample traits (using eigengene network methodology), 
#and for calculating module membership measures.)

#Weighted Gene Co-expression Network Analysis (WGCNA) allows one to find clusters (modules) of genes, 
#where each module has nodes (individual genes) whose expression profiles are closely interconnected. 
#Using modules, one can summarize the node profiles of a given module by looking at a hub gene, 
#which is a representative gene of the module since it is centrally located and highly interconnected with other genes in the module. 
#By comparing and relating modules instead of individual nodes to a sample trait, multiple testing problem is reduced. 
#One can also annotate how close a node is to identified modules to define its module membership to each module. 
#A hub gene would have a very high module membership to the module that it is the hub gene for. 
#Gene significance (GS) is a measure of how biological significant a given gene is to a sample trait, 
#where the higher the absolute value of GS is, the more biologically significant this gene is to the sample trait. 
#One can use GS to quantify the correlation strength of a gene expression profile to a sample trait. By finding GS, 
#one can incorporate external information into the co-expression network. 
#Module significance is the average absolute gene significance measure for all genes in a module. 
#So if there is a high number of genes with high GS in a module in relation to a sample trait, 
#then higher that module's module significance would be in relation to that sample trait. 
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559#citeas
#1. Do module-trait correlation to find trait of interest.
#2. Do barplot of trait-based mean gene significance across modules to find the module that has the highest gene significance to the trait. 
#3. A scatterplot of gene significance vs. module membership in the most significant module 
#will show which gene has the highest biological significance to a sample trait (GS) while also being the closest to its identified module (MM). 
#4. Do a hierarchical clustering dendrogram of module eigengenes and the sample trait of interest. 
#5. Heatmap plot of the adjacencies in the eigengene network including the trait y. 
#Each row and column in the heatmap corresponds to one module eigengene (labeled by color) or the trait (labeled by y). 
#In the heatmap, green color represents low adjacency (negative correlation), while red represents high adjacency (positive correlation).

#install stuff first
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute"))
BiocManager::install("WGCNA")
install.packages("BiocManager")
library(WGCNA)
library(flashClust)
library(cluster)
library(DOSE)
library(gplots)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#EV_Cell####
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/##

#let's get started
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

#
#Read in the female liver data set
#femData = read.csv("LiverFemale3600.csv");
# Read in the male liver data set
#maleData = read.csv("LiverMale3600.csv");
# Take a quick look at what is in the data sets (caution, longish output):
#dim(femData)
#names(femData)
#dim(maleData)
#names(maleData)

#1.Simulation of expression and trait data: PDF document, R script
#2.Loading of expression data, an alternative to data simulation, provided to illustrate data loading of real data: PDF document, R script
View(EV_vs_Cell_norm_counts_Gene_ID)
View(t_EV_vs_Cell_norm_counts_Gene_ID)

# the data needs the rows correspond to samples and the columns correspond to genes
rna_EV_Cell_norm_counts_wgcna<-t_EV_vs_Cell_norm_counts_Gene_ID
View(rna_EV_Cell_norm_counts_wgcna)


#3.Basic data pre-processing illustrates rudimentary techniques for handling missing data and removing outliers: PDF document, R script
#In this section we illustrate basic data cleaning and pre-processing steps for expression data.

#3.a Identification of outlying samples
#We start by determining the mean expression per array and the number of missing values per array:


meanExpressionByArray=apply(rna_EV_Cell_norm_counts_wgcna,1,mean, na.rm=T)

#for NumberMissingByArray,
#arrays with excessive numbers of missing data should be removed
#However, since we do not have any missing values but 0's, we will actually 
#assign NumberMissingByArray as NOT NumberMissingByArray=apply( is.na(data.frame(rna_EV_Cell_norm_counts_wgcna)),1, sum)
#but NumberMissingByArray=apply( data.frame(rna_EV_Cell_norm_counts_wgcna)==0,1, sum)

#NumberMissingByArray=apply( is.na(data.frame(rna_EV_Cell_norm_counts_wgcna)),1, sum)
NumberMissingByArray=apply( data.frame(rna_EV_Cell_norm_counts_wgcna)==0,1, sum)
View(NumberMissingByArray)
#then, arrays with excessive numbers of missing data should be removed
#While there are about 14000 to 23000 genes missing per sample/array out of 40000 possible genes. 
#we are not throwing away any sample since the missing numbers are similar. 
# For example, Keep only arrays containing less than 500 missing entries
#KeepArray= NumberMissingByArray<500
#table(KeepArray)


#datExpr=datExpr[KeepArray,]
#y=y[KeepArray]
#ArrayName[KeepArray]

#A simple way to examine the mean expression per array is to use
sizeGrWindow(9, 5)
barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:36), cex.names = 0.7)
#you end up with same mean expression across the 36 samples because they were normalized
#No arrays in the plot seem to have an outlying mean expression value. 
#The numbers of missing entries in each array are:
NumberMissingByArray
#While there are about 14000 to 23000 genes missing per sample/array out of 40000 possible genes. 
#we are not throwing away any sample since the missing numbers are similar. 

#3.b Handling missing data and zero variance in probe profiles
#Here we count the number of missing samples in each probe profile, and remove probes with extensive numbers of missing samples.
#However, since we do not have any missing values due to imputation, we will actually do
#no.present_rna_EV_Cell_norm_counts_wgcna=as.vector(apply(!(as.matrix(rna_EV_Cell_norm_counts_wgcna)==0),2, sum) )
#instead of no.present_rna_EV_Cell_norm_counts_wgcna=as.vector(apply(!is.na(as.matrix(rna_EV_Cell_norm_counts_wgcna)),2, sum) )

#In addition, we remove probes that do not vary at all.
NumberMissingByGene =apply( is.na(data.frame(rna_EV_Cell_norm_counts_wgcna)),2, sum)
# One could do a barplot(NumberMissingByGene), but the barplot is empty in this case.
# It may be better to look at the numbers of missing samples using the summary method:
summary(NumberMissingByGene)
# Calculate the variances of the probes and the number of present entries
variance_rna_EV_Cell_norm_counts_wgcna=as.vector(apply(as.matrix(rna_EV_Cell_norm_counts_wgcna),2,var, na.rm=T))
no.present_rna_EV_Cell_norm_counts_wgcna=as.vector(apply(!(as.matrix(rna_EV_Cell_norm_counts_wgcna)==0),2, sum) )

# Another way of summarizing the number of pressent entries
table(no.present_rna_EV_Cell_norm_counts_wgcna)
# Keep only genes whose variance is non-zero and have at least 27 present entries (ev plus cell equaling 27)
KeepGenes= variance_rna_EV_Cell_norm_counts_wgcna>0 & no.present_rna_EV_Cell_norm_counts_wgcna>=27
table(KeepGenes)
#> table(KeepGenes)
#KeepGenes
#FALSE  TRUE 
#20437 19768 

#So, 20437 proteins have 0 value in 9 or more samples, so we will exclude these. 
rna_EV_Cell_norm_counts_wgcna=rna_EV_Cell_norm_counts_wgcna[, KeepGenes]

#In this case, since the data is simulated without missing data or zero-variance probes, all probes are retained.
View(rna_EV_Cell_norm_counts_wgcna)
#3.c Rudimentary detection of outlier samples
#We use hierarchical clustering with the Euclidean distance to determine whether there are array (sample) outliers:
datExpr=rna_EV_Cell_norm_counts_wgcna
#let's change the rows/samples so that it is similar to proteomic data i have done for a previous paper.
#where it is ev first then cell, ev 3,4,5 female then ev 3,4,5 male, etc
datExpr<-datExpr[c(1:3,7:9,13:15,4:6,10:12,16:18,19:21,25:27,31:33,22:24,28:30,34:36),]

View(datExpr)
sizeGrWindow(9, 5)
plotClusterTreeSamples(datExpr)
#worked. so there is no value in x-axis, but the height/y-axis distance shows how close samples are to each other. 
#and samples in different branches, the distance they have to travel up the branch and down to the other branch to meet each other is 
#how far they are from each other. 
write.csv(rna_EV_Cell_norm_counts_wgcna,file="rna_EV_Cell_norm_counts_wgcna.csv")

#another way of doing Hierarchical clustering: 
dist <- dist(rna_EV_Cell_norm_counts_wgcna, diag=TRUE)
# Hierarchical Clustering with hclust
hc <- hclust(dist)
# Plot the result
plot(hc)
#worked in plotting the cluster dendrogram, but cannot know whether there are array/sample outliers. not necessary for my samples.
getwd()

#4.Standard gene screening illustrates gene selection based on Pearson correlation and shows that the results are not satisfactory: PDF document, R script
#not doing this one since we dont have clinical trait to match the genes



#5.Construction of a weighted gene co-expression network and network modules illustrated step-by-step; includes a discussion of alternate clustering techniques: PDF document, R script
# Load additional necessary packages
library(cluster)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#In this section we provide a step-by-step overview of gene network construction and module detection.

#5.a Defining a weighted gene co-expression network
# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1=abs(cor(datExpr,use="p"))^6
# When you have relatively few genes (<5000) use the following code. Use this if we have less than 4000 genes
#k is for connectivity
#k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
datE<-datExpr
k=softConnectivity(datE,power=6)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
#worked
#The left panel shows a histogram of network connectivities. The right panel shows a log-log plot of the same histogram.
#The approximate straight line relationship (with R2 value) shows approximate scale free topology. 
#In most applications we find that scale free topology is at least approximately
#satisfied when a high power is chosen for defining the adjacency matrix. 
#We should point out that is not necessary that a network satisfies scale free topology; 
#scale free topology may not be satisfied if the data are comprised of 
#globally very distinct groups of samples (e.g. different tissues types). 
#Poor fit to scale free topology may also indicate the presence of array outliers.
#My guess is that because ev and cell samples are so distinct from each other.

#5.c.2 Use of topologial overlap to define dissimilarity
#Adjacency can be used to define a separate measure of similarity, the Topological Overlap Matrix(TOM) [2, 1]:
dissTOM=TOMdist(ADJ1)
collectGarbage()

#let's get gene names from our data
datExpr
SubGeneNames=colnames(datExpr)

#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#choose softpower
#If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers 
#(less than 15 for unsigned or signed hybrid networks, and less than 30 for signed networks) 
#and the mean connectivity remains relatively high (in the hundreds or above), 
#chances are that the data exhibit a strong driver that makes a subset of the samples globally different from the rest. 
#The difference causes high correlation among large groups of genes 
#which invalidates the assumption of the scale-free topology approximation.

#If the lack of scale-free topology fit turns out to be caused by an interesting biological variable 
#that one does not want to remove (i.e., adjust the data for), 
#the appropriate soft-thresholding power can be chosen based on the number of samples as in the table below. 
#This table has been updated in December 2017 to make the resulting networks conservative.
#Number of samples	Unsigned and signed hybrid networks	Signed networks
#Less than 20	                  9	                            18
#20-30	                        8	                            16
#30-40	                        7	                            14
#more than 40	                  6	                            12

powers = c(c(1:10), seq(from = 12, to=20, by=2));
softPower=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "signed")
#
#Plot the results
sft<-softPower

#to get 2 separate plots in one window:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#so, looking at the mean connectivity, you usually want to choose a soft threshold (power) that is over 0.8 for R^2,
#which would be 14 for 36 EV and cell samples, so lets try that.
#for EV, it is 6 and for Cell, it is 9. The reason 36 samples of EV and Cell is 14 is due to the fact that 
#EV and Cell sample groups are so different from each other. 
#so use softPower 14 for the 36 samples, and use softPower 9 for EV and Cell so that we can compare EV with Cell. 
#for more info: https://bioinformatics.stackexchange.com/questions/11335/wgcna-problem-with-selecting-soft-threshold
softPower = 14

#let's do a adjacency matrix (kind of like a correlation matrix)
adj= adjacency(datExpr,type = "signed", power = softPower);

#turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
#trying to calculate how connected genes are to each other.
#we are using signed for network and TOM type instead of unsigned.
#Peter Langfelder recommends using signed network.
#https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/
#How should pairs of nodes with strong negative correlations be treated in a correlation network analysis? 
#One option is to consider them connected, just as if the correlation were positive. 
#A network constructed in this way is an unsigned network, because the sign of the correlation does not matter. 
#On the other hand, strongly negatively correlated nodes can also be considered unconnected. This leads to a signed network, 
#so called because the sign of a strong correlation value makes all the difference 
#between the pair of nodes being strongly connected or not connected at all. 
#To avoid any confusion, I want to emphasize that the resulting adjacency matrix 
#(the matrix that contains the connection strengths between nodes) is always non-negative.
#First, more often than not, direction does matter: 
#it is important to know where node profiles go up and where they go down, 
#and mixing negatively correlated nodes together necessarily mixes the two directions together. 
#Second, negatively correlated nodes often belong to different categories. 
#For example, in gene expression data, negatively correlated genes tend to come from biologically very different categories. 
View(datExpr)
# Convert all variable types to numeric
datExpr <- as.data.frame(apply(datExpr, 2, as.numeric))  

TOM=TOMsimilarityFromExpr(datExpr,networkType = "signed", TOMType = "signed", power = softPower);
#since correlation matrix will not have names, we are giving subgenenames to column and row names. 
colnames(TOM) =rownames(TOM) =SubGeneNames
dissTOM=1-TOM
#
#hierarchical clustering of the genes based on the TOM dissimilarity measure
library(flashClust)
geneTree = flashClust(as.dist(dissTOM),method="average")

#plot the resulting clustering tree (dendrogram)
par(mfrow = c(1,1))
plot(geneTree, xlab="", sub="",cex=0.3,main=paste0("genes"))
#from that, we see how the genes are clustered together. 
#we will cut the tree at a horizontal line, let's say 0.8 height, then whatever is cut below will form a module.
#dynamiccuttree will do this automatically for you.

# Set the minimum module size
minModuleSize = 20;

# Module identification using dynamic tree cut
#important to identify modules
dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
#static
#colorStaticADJ=as.character(cutreeStaticColor(geneTree, cutHeight=.975, minSize=20))
# Plot the dendrogram with module colors

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
#tells how many genes i have in each module and how many modules there are. 
table(dynamicMods)
#table(colorStaticADJ)
##Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#staticColors = labels2colors(colorStaticADJ)
#table(staticColors)
#
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
#so grey color is leftover genes that dont really belong to any module
#the other modules are important, which we will extract and study more later. 
#

#dynamicColors=staticColors
#below will give you module tables with gene names in each module
module_colors= setdiff(unique(dynamicColors), "grey")

for (color in module_colors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_EV_Cell",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}
#Dae, run colorh1= dynamicColors
colorh1= dynamicColors
#Dae, to get hub genes, 
hubs    = chooseTopHubInEachModule(datExpr, colorh1)
hubs

#for Gene_ID and moduleColor, that set of codes were run few hundred lines below with geneInfo0


###TOM heatmap plot###
library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
#dynamicColors=staticColors
#set the diagonal of the dissimilarity to NA 
diag(dissTOM) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to a power  to bring out the module structure
#8.b Topological overlap matrix plot for visualizing the network
#We now create a so-called TOM plot, a heatmap plot depicting the topological overlap matrix 
#supplemented by hierarchical clustering dendrograms and the module colors.

png("rna_WGCNA_EV_Cell_TOM Heatmap Plot, Module RNA_sft_pwr_14.png", width = 2500, height = 2500, res = 300)
TOMplot(dissTOM^4, geneTree, as.character(dynamicColors), main = "EV and Cell TOM Heatmap Plot, Module RNAs", col=myheatcol )
dev.off()


#let's do a multidimensional scaling plot
#Multidimensional scaling (MDS) is a multivariate data analysis approach 
#that is used to visualize the similarity/dissimilarity between samples by plotting points in two dimensional plots.
cmd1=cmdscale(as.dist(dissTOM),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(colorh1), main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")



####modules to external traits relationship###
#here, we want to see which modules relate closely to which traits
#for my samples, traits will be pregnancy/set, sex, and ethanol.
#it can also be ev or cell, if i have all the samples together.
#Module-Trait relationships. Color scale (red-green) represents the strength of the correlation between the module and the trait. 
#Each box gives a correlation value (R^2) followed by p-value (in parenthesis). 
#the color scale is the strength of the correlation between the module and the trait, 
#where for example, MEblack and Female being very red means highly correlated.

#let's load the excel csv file i have created for traits. 
#i am adding "row.names=1" function, which means make the first column (sample names) into row names 
#datTraits = read.csv("Traits_wgcna1.csv",row.names = 1)
#write.csv(datTraits,file="Traitswgcna.csv")
#form a data frame analogous to expression data that will hold the clinical traits.
#let's match the rownames of datExpr and datTraits, after you check that the samples are lined up correctly by rows
rownames(datExpr)<-rownames(datTraits)

table(rownames(datTraits)==rownames(datExpr)) 
#should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)

##
getwd()

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


#3 Relating modules to external clinical traits
#3.a Quantifying module-trait associations
#In this analysis we would like to identify modules that are significantly associated with the measured clinical traits.
#Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external
#traits and look for the most significant associations:

#instead of moduleColors, i have dynamicColors from:
#dynamicColorsEV = labels2colors(dynamicModsEV)
View(dynamicColors)
View(colorh1)

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
#moduleEigengenes() simply calculate the 1st Principal Component (PC) i.e., module eigengene (ME), of each module.
MEs0 = moduleEigengenes(datExpr, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#let's use a graphical representation to read the table, a color-coded table,
#where we color code each association of module eigengenes to the sample trait by the correlation value:
sizeGrWindow(16,10)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))

#View(datTraits)
#colnames(datTraits)[1]<-"EV_Cell"
#colnames(datTraits)[2]<-"Pregnancy"
#colnames(datTraits)[4]<-"Ethanol"
#colnames(datTraits)[1]<-"EV_Cell"
#looking at module-trait relationships, there is a significant strong positive correlation of MEbrown, blue, turquoise modules to Cell trait,
#while rest are for EV trait. 


#Do barplot of trait-based mean gene significance across modules to find the module that has the highest gene significance to the trait. 
#III. Using simulated data to evaluate different module detection methods
#and gene screening approaches
#6. Relating modules and module eigengenes to external data
#6.a Representing modules by eigengenes and relating eigengenes to one another
#To get a sense of how related the modules are one can summarize each module by its eigengene (first principal component).


datME=moduleEigengenes(datExpr,colorh1)$eigengenes
signif(cor(datME, use="p"), 2)
#The result is

#
#MEblack MEblue MEbrown MEcyan MEdarkred MEgreen MEgreenyellow MEgrey60 MElightcyan MElightgreen MElightyellow
#MEblack          1.000 -0.710   -0.22  0.850    0.3100   0.890         0.990    0.350        0.69         0.82       -0.5000
#MEblue          -0.710  1.000    0.62 -0.680   -0.4800  -0.670        -0.720   -0.630       -0.72        -0.73       -0.0930
#MEbrown         -0.220  0.620    1.00 -0.240   -0.2200  -0.210        -0.240   -0.260       -0.30        -0.28       -0.4600
#MEcyan           0.850 -0.680   -0.24  1.000    0.4900   0.850         0.890    0.240        0.82         0.82       -0.3100
#MEdarkred        0.310 -0.480   -0.22  0.490    1.0000   0.450         0.260    0.330        0.85         0.78        0.0015
#MEgreen          0.890 -0.670   -0.21  0.850    0.4500   1.000         0.870    0.360        0.69         0.84       -0.3900
#MEgreenyellow    0.990 -0.720   -0.24  0.890    0.2600   0.870         1.000    0.360        0.69         0.78       -0.4600
#MEgrey60         0.350 -0.630   -0.26  0.240    0.3300   0.360         0.360    1.000        0.46         0.43        0.0180
#MElightcyan      0.690 -0.720   -0.30  0.820    0.8500   0.690         0.690    0.460        1.00         0.92       -0.1900
#MElightgreen     0.820 -0.730   -0.28  0.820    0.7800   0.840         0.780    0.430        0.92         1.00       -0.3100
#MElightyellow   -0.500 -0.093   -0.46 -0.310    0.0015  -0.390        -0.460    0.018       -0.19        -0.31        1.0000
#MEmagenta        0.540 -0.790   -0.32  0.380    0.4100   0.460         0.520    0.860        0.59         0.59       -0.0780
#MEmidnightblue   0.510 -0.620   -0.25  0.640    0.9500   0.670         0.470    0.420        0.89         0.89       -0.0910
#MEpink           0.630 -0.810   -0.34  0.600    0.7700   0.590         0.600    0.690        0.88         0.86       -0.1100
#MEpurple         0.820 -0.690   -0.26  0.780    0.2500   0.540         0.860    0.380        0.70         0.64       -0.3300
#MEred            0.770 -0.450   -0.11  0.700    0.1800   0.930         0.750    0.160        0.40         0.61       -0.4200
#MEroyalblue      0.120 -0.240   -0.10  0.057    0.0650   0.130         0.130    0.590        0.15         0.13        0.0580
#MEsalmon         0.830 -0.700   -0.24  0.540    0.3200   0.600         0.780    0.520        0.64         0.71       -0.4400
#MEtan           -0.600  0.044   -0.26 -0.410   -0.0750  -0.480        -0.560   -0.037       -0.30        -0.41        0.9700
#MEturquoise     -0.740  0.850    0.12 -0.690   -0.4700  -0.690        -0.740   -0.620       -0.72        -0.74        0.1800
#MEyellow        -0.016 -0.570   -0.32  0.120    0.2400   0.049         0.028    0.430        0.23         0.14        0.7200
#MEmagenta MEmidnightblue MEpink MEpurple  MEred MEroyalblue MEsalmon  MEtan MEturquoise MEyellow
#MEblack            0.540          0.510   0.63     0.82  0.770       0.120    0.830 -0.600       -0.74   -0.016
#MEblue            -0.790         -0.620  -0.81    -0.69 -0.450      -0.240   -0.700  0.044        0.85   -0.570
#MEbrown           -0.320         -0.250  -0.34    -0.26 -0.110      -0.100   -0.240 -0.260        0.12   -0.320
#MEcyan             0.380          0.640   0.60     0.78  0.700       0.057    0.540 -0.410       -0.69    0.120
#MEdarkred          0.410          0.950   0.77     0.25  0.180       0.065    0.320 -0.075       -0.47    0.240
#MEgreen            0.460          0.670   0.59     0.54  0.930       0.130    0.600 -0.480       -0.69    0.049
#MEgreenyellow      0.520          0.470   0.60     0.86  0.750       0.130    0.780 -0.560       -0.74    0.028
#MEgrey60           0.860          0.420   0.69     0.38  0.160       0.590    0.520 -0.037       -0.62    0.430
#MElightcyan        0.590          0.890   0.88     0.70  0.400       0.150    0.640 -0.300       -0.72    0.230
#MElightgreen       0.590          0.890   0.86     0.64  0.610       0.130    0.710 -0.410       -0.74    0.140
#MElightyellow     -0.078         -0.091  -0.11    -0.33 -0.420       0.058   -0.440  0.970        0.18    0.720
#MEmagenta          1.000          0.530   0.87     0.58  0.200       0.270    0.750 -0.160       -0.79    0.430
#MEmidnightblue     0.530          1.000   0.84     0.35  0.430       0.130    0.450 -0.180       -0.62    0.260
#MEpink             0.870          0.840   1.00     0.65  0.270       0.200    0.770 -0.210       -0.81    0.380
#MEpurple           0.580          0.350   0.65     1.00  0.320       0.120    0.790 -0.420       -0.71    0.130
#MEred              0.200          0.430   0.27     0.32  1.000       0.079    0.390 -0.490       -0.47   -0.110
#MEroyalblue        0.270          0.130   0.20     0.12  0.079       1.000    0.170  0.047       -0.24    0.250
#MEsalmon           0.750          0.450   0.77     0.79  0.390       0.170    1.000 -0.530       -0.72    0.037
#MEtan             -0.160         -0.180  -0.21    -0.42 -0.490       0.047   -0.530  1.000        0.22    0.710
#MEturquoise       -0.790         -0.620  -0.81    -0.71 -0.470      -0.240   -0.720  0.220        1.00   -0.520
#MEyellow           0.430          0.260   0.38     0.13 -0.110       0.250    0.037  0.710       -0.52    1.000


#We define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation
#between the module eigengenes, and use it to cluster the eigengene:
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
sizeGrWindow(8,9)
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes (first PC)")

#6.a.1 Pairwise scatter plots of the samples (arrays) along the module eigengenes
#We create a pairwise scatter plots of the samples (arrays) along the module eigengenes:
#define the microarray sample trait
y=datTraits$EV_Cell


sizeGrWindow(8,9)
plotMEpairs(datME,y=y)

#6.b Relating observed module eigengenes to simulated true module eigengenes
#We now relate the observed module eigengenes to the simulated true module eigengenes:
#attach(ModuleEigengeneNetwork1)
#ModuleEigengeneNetwork1=data.frame(y,MEturquoise,MEblue,MEbrown,MEyellow)

#6.e.2 Measure of module significance as average gene significance
#One can also define a measure of module significance as the average gene significance of all genes in the module. We
#use the absolute value for defining a correlation based gene significance measure.
GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, dissimME=(1-t(cor(datME, method="p")))/2
                          , mean, na.rm=T)
ModuleSignificance
#here is the result:
#black           blue          brown           cyan      darkgreen       darkgrey darkolivegreen 
#0.3745184      0.9010115      0.3167285      0.4541212      0.4287165      0.4695732      0.4124254 
#darkorange        darkred  darkturquoise          green    greenyellow           grey         grey60 
#0.4053298      0.4838315      0.4937343      0.3161906      0.5248919      0.1188808      0.4036407 
#lightcyan     lightgreen    lightyellow        magenta   midnightblue         orange  paleturquoise 
#0.4923496      0.5595521      0.5760983      0.2962100      0.5044114      0.5484404      0.3906726 
#pink         purple            red      royalblue    saddlebrown         salmon        skyblue 
#0.7438960      0.5791215      0.5649846      0.6061930      0.5883377      0.4087688      0.4589132 
#steelblue            tan      turquoise         violet          white         yellow 
#0.3714220      0.5277154      0.7563999      0.2351134      0.4663825      0.4679870 

#To plot module significance, one can use
sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,colorh1)

#The advantage of this approach is that it can be used for any gene significance measure. 
#A gene significance measure could be defined without reference to a sample trait. 
#For example, it could indicate pathway membership (1 or 0) or gene essentiality (1 or 0), etc.
#The next logical step in an analysis of empirical data would be to carry out a functional enrichment analysis of
#each module, for example using the software EASE (David) at http://david.abcc.ncifcrf.gov/summary.jsp.
#We refer the reader to Tutorial I for an example of a functional enrichment analysis. Before we end, we again save
#calculated results for use in subsequent sections.
collectGarbage()
save.image("Transcriptomic_WGCNA_v4.RData")


#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf

#3.b Gene relationship to trait and important modules: Gene Significance and Module Membership
#We quantify associations of individual genes with our trait of interest (EV_Cell) by defining Gene Significance GS as
#(the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative
#measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 
#This allows us to quantify the similarity of all genes on the array to every module.
# Define variable trait containing the EV_Cell column of datTrait
EV_Cell = as.data.frame(datTraits$EV_Cell);
names(EV_Cell) = "EV_Cell"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, EV_Cell, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(EV_Cell), sep="");
names(GSPvalue) = paste("p.GS.", names(EV_Cell), sep="");


#3.c Intramodular analysis: identifying genes with high GS and MM
#Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module
#membership in interesting modules. As an example, we look at the blue module that has the highest association
#with EV_Cell We plot a scatterplot of Gene Significance vs. Module Membership in the blue module:
#module blue is highest gene significance for Cell,_____ for EV.
# names (colors) of the modules
library(ggplot2)

module = "blue"
column = match(module, modNames);
moduleGenes = dynamicColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for EV_Cell location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "turquoise"
column = match(module, modNames);
moduleGenes = dynamicColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for EV_Cell location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "magenta"
column = match(module, modNames);
moduleGenes = dynamicColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for EV_Cell location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module)


module = "yellow"
column = match(module, modNames);
moduleGenes = dynamicColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for EV_Cell location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



#I. Network analysis of liver expression data in female mice
#5. Network visualization using WGCNA functions
#5.b Visualizing the network of eigengenes
#It is often interesting to study the relationships among the found modules. 
#One can use the eigengenes as representative profiles and quantify module similarity by eigengene correlation. 
#The package contains a convenient function plotEigengeneNetworks that generates a summary plot of the eigengene network. 
#It is usually informative to add a clinical trait (or multiple traits) to the eigengenes 
#to see how the traits fit into the eigengene network:
#Previously: 
#MEs0 = moduleEigengenes(datExpr, dynamicColors)$eigengenes
#MEs = orderMEs(MEs0)
## Isolate weight from the clinical traits
#EV_Cell = as.data.frame(datTraits$EV_Cell);
#names(EV_Cell) = "EV_Cell"
# Add the EV_Cell to existing module eigengenes
MET = orderMEs(cbind(MEs, EV_Cell))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
#The function produces a dendrogram of the eigengenes and trait(s), and a heatmap of their relationships. 
#To split the dendrogram and heatmap plots, we can use the following code
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

collectGarbage()


#3.d Summary output of network analysis results
#We have found modules with high association with our trait of interest, and have identified their central players by
#the Module Membership measure. We now merge this statistical information with gene annotation and write out a
#file that summarizes the most important results and can be inspected in standard spreadsheet software such as MS
#Excel or Open Office Calc. Our expression data are only annotated by probe ID names: the command
colnames(datExpr)
#will return all Gene IDs included in the analysis. Similarly,
colnames(datExpr)[dynamicColors=="blue"]
#will return probe IDs belonging to the blue module. To facilitate interpretation of the results, 
#we will connect Uniprot IDs to gene names and universally recognized identification numbers (Entrez codes).
View(datExpr)
#will need to convert to Enterz ID
#https://yulab-smu.github.io/clusterProfiler-book/chapter14.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
#install db for mouse IDs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DOSE")
library(DOSE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ReactomePA")
library(ReactomePA)

getwd()
#create vector of IDs
View(datExpr)
tdatExpr<-as.data.frame(t(datExpr))
View(tdatExpr)
library(data.table)
setDT(tdatExpr, keep.rownames = "Gene_ID")[]

ID <- tdatExpr$Gene_ID
eg <- bitr(ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg)
#rename column in eg so matching title
colnames(eg)[1] <- "Gene_ID"
#merge based on Gene_ID
tdatExpr_Entrez <- merge(eg, tdatExpr, all.y = TRUE, by = "Gene_ID")
View(tdatExpr_Entrez)
View(geneTraitSignificance)
View(GSPvalue)

#because 4.86% of input gene IDs are fail to map, we lost some proteins, going from 2644 to 2615 (old one was 2238 to 2124). 
#and since geneTraitSignificance and GSPvalue both have 2238, we need to merge them to a new d.f. with 2644 rows.



#We now create a data frame holding the following information for all probes: probe ID, gene symbol, Locus Link ID
#(Entrez code), module color, gene significance for EV_Cell, and module membership and p-values in all modules. The
#modules will be ordered by their significance for EV_Cell, with the most significant ones to the left.
# Create the starting data frame


#let's start#
geneInfo0<-data.frame(Gene_ID = tdatExpr$Gene_ID,
                      ENTREZID = tdatExpr_Entrez$ENTREZID,
                      moduleColor = dynamicColors)

#previously done above:
EV_Cell = as.data.frame(datTraits$EV_Cell);
names(EV_Cell) = "EV_Cell"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, EV_Cell, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(EV_Cell), sep="");
names(GSPvalue) = paste("p.GS.", names(EV_Cell), sep="");
#up to here.


# Order modules by their significance for sample trait EV_Cell location
modOrder = order(-abs(cor(MEs, EV_Cell, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

#then let's get each gene's correlation and significance to the trait of choice (Location: EV vs cell)
geneTraitSignificance1<-geneTraitSignificance
GSPvalue1<-GSPvalue
setDT(geneTraitSignificance1, keep.rownames = "Gene_ID")[]
View(geneTraitSignificance1)
setDT(GSPvalue1, keep.rownames = "Gene_ID")[]

#then merge. for whatever reason, when i merged first then did data.frame() to merge different dfs, 
#i would end up with rows shifting during merging and ending with wrong values in rows. it was "wrowng"
#thus, we do data.frame() first, then merge() by Gene_ID
geneInfo0<-merge(geneInfo0,geneTraitSignificance1, by.y = "Gene_ID")
geneInfo0<-merge(geneInfo0,GSPvalue1, by.y = "Gene_ID")
#once i merge, check to see if i have different number of rows
#let's see how many unique rows i have in Protein column:
length(unique(geneInfo0$Gene_ID))
#19768 
#if i have duplicates,
#remove them. 
geneInfo0<-geneInfo0[!duplicated(geneInfo0$Gene_ID), ]
#geneInfo0<-geneInfo0[,c(1,2,30,31,3:29)]
View(geneInfo0)
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.EV_Cell));
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
getwd()

write.csv(geneInfo, file = "rna_WGCNA_EV_Cell_geneInfo_sft_pwr_14.csv")

View(geneInfo)

#let's pick columns with Gene_ID, ENTREZID, and moduleColor
rna_WGCNA_EV_Cell_Gene_ID_moduleColor_sft_pwr_14<-
  geneInfo[,c("Gene_ID","ENTREZID","moduleColor")]
write.csv(rna_WGCNA_EV_Cell_Gene_ID_moduleColor_sft_pwr_14, file = "rna_WGCNA_EV_Cell_Gene_ID_moduleColor_sft_pwr_14.csv")



collectGarbage()

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute"))
BiocManager::install("WGCNA")
install.packages("BiocManager")
library(WGCNA)
library(flashClust)
library(cluster)
library(DOSE)
library(gplots)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#4 Interfacing network analysis with other data such as functional annotation and gene ontology
#Our previous analysis has identified several modules (labeled brown, red, and salmon) that are highly associated with weight. 
#To facilitate a biological interpretation, we would like to know the gene ontologies of the genes in the modules, 
#whether they are significantly enriched in certain functional categories etc.

#4.a Output gene lists for use with online software and services

#One option is to simply export a list of gene identifiers that can be used as input for several popular gene ontology
#and functional enrichment analysis suites such as David or AmiGO. 
#For example, we write out the LocusLinkID (entrez) codes for the brown module into a file:

View(datExpr)
# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExpr)
#probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = geneInfo$ENTREZID;
# $ Choose interesting modules
intModules = c("blue", "brown", "yellow", "green", "turquoise")
for (module in intModules)
{
  # Select module probes
  modGenes = (dynamicColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
#probes = names(datExpr)
#probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
#allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
#intModules = c("brown", "red", "salmon")
#for (module in intModules)
#{
# Select module probes
#  modGenes = (moduleColors==module)
# Get their entrez ID codes
#  modLLIDs = allLLIDs[modGenes];
# Write them into a file
#  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
#  write.table(as.data.frame(modLLIDs), file = fileName,
#              row.names = FALSE, col.names = FALSE)
#}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)
View(allLLIDs)


GOenr = GOenrichmentAnalysis(dynamicColors, allLLIDs, organism = "mouse", nBestP = 10);

tab = GOenr$bestPTerms[[4]]$enrichment

View(tab)
names(tab)

write.table(tab, file = "rna_WGCNA_EV_Cell_GOEnrichmentTable_sft_pwer_14.csv", sep = ",", quote = TRUE, row.names = FALSE)



keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab
View(screenTab)
View(tab)
write.table(screenTab, file = "rna_WGCNA_EV_Cell_GOEnrichmentTable_shortened_sft_pwr_14.csv", sep = ",", quote = TRUE, row.names = FALSE)

save.image("Transcriptomic_WGCNA_v4.RData")
getwd()




###let's use merged tree cut instead of dynamic tree cut####
#2.b.5 Merging of modules whose expression profiles are very similar
#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge
#such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we
#calculate their eigengenes and cluster them on their correlation:
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
sizeGrWindow(8,9)
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes (first PC)")

#We choose a height cut of 0.2, corresponding to correlation of 0.8, to merge (see Fig. 4):
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
#To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged
#module colors underneath (Figure 5).
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
table(dynamicColors)
table(mergedColors)
png("rna_WGCNA_EV_Cell_TOM Heatmap Plot_mergedColors_h0.2, Module RNA_sft_pwr_14.png", width = 2500, height = 2500, res = 300)
TOMplot(dissTOM^4, geneTree, as.character(mergedColors), main = "EV and Cell TOM Heatmap Plot, Module RNAs", col=myheatcol )
dev.off()

#dynamicColors=staticColors
#below will give you module tables with gene names in each module
module_colors= setdiff(unique(mergedColors), "grey")

for (color in module_colors){
  module=SubGeneNames[which(mergedColors==color)]
  write.table(module, paste("module_EV_Cell",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}




mergedcolorh1= mergedColors
#Dae, to get hub genes, 
mergedhubs    = chooseTopHubInEachModule(datExpr, mergedcolorh1)
mergedhubs

table(mergedcolorh1)
mergedModuleSignificance


cmd1=cmdscale(as.dist(dissTOM),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(mergedcolorh1), main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

# Recalculate MEs with color labels
#moduleEigengenes() simply calculate the 1st Principal Component (PC) i.e., module eigengene (ME), of each module.
mergedMEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
mergedMEs = orderMEs(mergedMEs0)
mergedmoduleTraitCor = cor(mergedMEs, datTraits, use = "p");
mergedmoduleTraitPvalue = corPvalueStudent(mergedmoduleTraitCor, nSamples);

#let's use a graphical representation to read the table, a color-coded table,
#where we color code each association of module eigengenes to the sample trait by the correlation value:
sizeGrWindow(16,10)
# Will display correlations and their p-values
textMatrix = paste(signif(mergedmoduleTraitCor, 2), "\n(",
                   signif(mergedmoduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(mergedmoduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = mergedmoduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))

#let's match the rownames of datExpr and datTraits, after you check that the samples are lined up correctly by rows
rownames(datExpr)<-rownames(datTraits)

table(rownames(datTraits)==rownames(datExpr)) 
#should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)
y=datTraits$EV_Cell

y
mergeddatME=moduleEigengenes(datExpr,mergedcolorh1)$eigengenes
signif(cor(mergeddatME, use="p"), 2)

#6.e.2 Measure of module significance as average gene significance
#One can also define a measure of module significance as the average gene significance of all genes in the module. We
#use the absolute value for defining a correlation based gene significance measure.
GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
mergedModuleSignificance=tapply(GeneSignificance, mergedcolorh1, mean, na.rm=T)


mergedModuleSignificance
#here is the result:
#        black          blue         brown    darkorange darkturquoise         green   greenyellow          grey     lightcyan        salmon 
#0.3745184     0.8153419     0.3167285     0.4698657     0.3880264     0.3359671     0.5328016     0.1188808     0.5365143     0.4087688 
#tan        violet        yellow 
#0.5277154     0.2351134     0.4679870 

#To plot module significance, one can use
sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,mergedcolorh1)

#We define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation
#between the module eigengenes, and use it to cluster the eigengene:
mergeddissimME=(1-t(cor(mergeddatME, method="p")))/2
mergedhclustdatME=hclust(as.dist(mergeddissimME), method="average" )
# Plot the eigengene dendrogram
sizeGrWindow(8,9)
par(mfrow=c(1,1))
plot(mergedhclustdatME, main="Clustering tree based of the module eigengenes (first PC)")


#The advantage of this approach is that it can be used for any gene significance measure. 
#A gene significance measure could be defined without reference to a sample trait. 
#For example, it could indicate pathway membership (1 or 0) or gene essentiality (1 or 0), etc.
#The next logical step in an analysis of empirical data would be to carry out a functional enrichment analysis of
#each module, for example using the software EASE (David) at http://david.abcc.ncifcrf.gov/summary.jsp.
#We refer the reader to Tutorial I for an example of a functional enrichment analysis. Before we end, we again save
#calculated results for use in subsequent sections.
collectGarbage()
save.image("Transcriptomic_WGCNA_v4.RData")


#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf

#3.b Gene relationship to trait and important modules: Gene Significance and Module Membership
#We quantify associations of individual genes with our trait of interest (EV_Cell) by defining Gene Significance GS as
#(the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative
#measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 
#This allows us to quantify the similarity of all genes on the array to every module.
# Define variable trait containing the EV_Cell column of datTrait
EV_Cell = as.data.frame(datTraits$EV_Cell);
names(EV_Cell) = "EV_Cell"
# names (colors) of the modules
modNames = substring(names(mergedMEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, mergedMEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, EV_Cell, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(EV_Cell), sep="");
names(GSPvalue) = paste("p.GS.", names(EV_Cell), sep="");


#3.c Intramodular analysis: identifying genes with high GS and MM
#Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module
#membership in interesting modules. As an example, we look at the blue module that has the highest association
#with EV_Cell We plot a scatterplot of Gene Significance vs. Module Membership in the blue module:
#module blue is highest gene significance for Cell,_____ for EV.
# names (colors) of the modules
library(ggplot2)

module = "blue"
column = match(module, modNames);
moduleGenes = dynamicColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for EV_Cell location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "brown"
column = match(module, modNames);
moduleGenes = dynamicColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for EV_Cell location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "magenta"
column = match(module, modNames);
moduleGenes = dynamicColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for EV_Cell location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module)


module = "lightcyan"
column = match(module, modNames);
moduleGenes = mergedColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for EV_Cell location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



#I. Network analysis of liver expression data in female mice
#5. Network visualization using WGCNA functions
#5.b Visualizing the network of eigengenes
#It is often interesting to study the relationships among the found modules. 
#One can use the eigengenes as representative profiles and quantify module similarity by eigengene correlation. 
#The package contains a convenient function plotEigengeneNetworks that generates a summary plot of the eigengene network. 
#It is usually informative to add a clinical trait (or multiple traits) to the eigengenes 
#to see how the traits fit into the eigengene network:
#Previously: 
#MEs0 = moduleEigengenes(datExpr, dynamicColors)$eigengenes
#MEs = orderMEs(MEs0)
## Isolate weight from the clinical traits
#EV_Cell = as.data.frame(datTraits$EV_Cell);
#names(EV_Cell) = "EV_Cell"
# Add the EV_Cell to existing module eigengenes
MET = orderMEs(cbind(mergedMEs, EV_Cell))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
#The function produces a dendrogram of the eigengenes and trait(s), and a heatmap of their relationships. 
#To split the dendrogram and heatmap plots, we can use the following code
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

collectGarbage()


#3.d Summary output of network analysis results
#We have found modules with high association with our trait of interest, and have identified their central players by
#the Module Membership measure. We now merge this statistical information with gene annotation and write out a
#file that summarizes the most important results and can be inspected in standard spreadsheet software such as MS
#Excel or Open Office Calc. Our expression data are only annotated by probe ID names: the command
colnames(datExpr)
#will return all Gene IDs included in the analysis. Similarly,
colnames(datExpr)[dynamicColors=="blue"]
#will return probe IDs belonging to the blue module. To facilitate interpretation of the results, 
#we will connect Uniprot IDs to gene names and universally recognized identification numbers (Entrez codes).

#will need to convert to Enterz ID
#https://yulab-smu.github.io/clusterProfiler-book/chapter14.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
#install db for mouse IDs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DOSE")
library(DOSE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ReactomePA")
library(ReactomePA)

getwd()
#create vector of IDs
View(datExpr)
tdatExpr<-as.data.frame(t(datExpr))
View(tdatExpr)
library(data.table)
setDT(tdatExpr, keep.rownames = "Gene_ID")[]

ID <- tdatExpr$Gene_ID
eg <- bitr(ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg)
#rename column in eg so matching title
colnames(eg)[1] <- "Gene_ID"
#merge based on Gene_ID
tdatExpr_Entrez <- merge(eg, tdatExpr, all.y = TRUE, by = "Gene_ID")
View(tdatExpr_Entrez)
View(geneTraitSignificance)
View(GSPvalue1)


#i put in the mergedColors as moduleColor here, because the moduleColor gets downshifted with geneInfo0 and does not
#match correctly if i do it later. 
#i found it out by catching that hub gene's module color did not match when looking at geneInfo0
#colnames(datExpr)[mergedColors=="lightcyan"]
#will return probe IDs belonging to the lightcyan module. 
#which(colnames(datExpr)[mergedColors=="lightcyan"]=="Ush2a")
#will tell you the location of Ush2a from the lightcyan module in datExpr. 
#If Ush2a is not in that named module, then it will return the value of integer(0)
#which(colnames(datExpr)[mergedColors=="blue"]=="Kpna1")

#let's start
geneInfo0<-data.frame(Gene_ID = tdatExpr$Gene_ID,
                      ENTREZID = tdatExpr_Entrez$ENTREZID,
                      moduleColor = mergedColors)
View(geneInfo0)

#We now create a data frame holding the following information for all probes: probe ID, gene symbol, Locus Link ID
#(Entrez code), module color, gene significance for EV_Cell, and module membership and p-values in all modules. The
#modules will be ordered by their significance for EV_Cell, with the most significant ones to the left.
# Create the starting data frame

#previously done above:
EV_Cell = as.data.frame(datTraits$EV_Cell);
names(EV_Cell) = "EV_Cell"
# names (colors) of the modules
modNames = substring(names(mergedMEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, mergedMEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, EV_Cell, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(EV_Cell), sep="");
names(GSPvalue) = paste("p.GS.", names(EV_Cell), sep="");
#up to here.


# Order modules by their significance for sample trait EV_Cell location
modOrder = order(-abs(cor(mergedMEs, EV_Cell, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

#then let's get each gene's correlation and significance to the trait of choice (Location: EV vs cell)
geneTraitSignificance1<-geneTraitSignificance
GSPvalue1<-GSPvalue
setDT(geneTraitSignificance1, keep.rownames = "Gene_ID")[]
View(geneTraitSignificance1)
setDT(GSPvalue1, keep.rownames = "Gene_ID")[]

#then merge. for whatever reason, when i merged first then did data.frame() to merge different dfs, 
#i would end up with rows shifting during merging and ending with wrong values in rows. it was "wrowng"
#thus, we do data.frame() first, then merge() by Gene_ID
geneInfo0<-merge(geneInfo0,geneTraitSignificance1, by.y = "Gene_ID")
geneInfo0<-merge(geneInfo0,GSPvalue1, by.y = "Gene_ID")
#once i merge, check to see if i have different number of rows
#let's see how many unique rows i have in Protein column:
length(unique(geneInfo0$Gene_ID))
#19768 
#if i have duplicates,
#remove them. 
geneInfo0<-geneInfo0[!duplicated(geneInfo0$Gene_ID), ]
geneInfo0<-geneInfo0[,c(1,2,30,31,3:29)]
View(geneInfo0)
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.EV_Cell));
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
getwd()

write.csv(geneInfo, file = "rna_WGCNA_EV_Cell_geneInfo_sft_pwr_14_h0.2.csv")

View(geneInfo)

#let's pick columns with Gene_ID, ENTREZID, and moduleColor
rna_WGCNA_EV_Cell_Gene_ID_moduleColor_sft_pwr_14_h0.2<-
  geneInfo[,c("Gene_ID","ENTREZID","moduleColor")]
write.csv(rna_WGCNA_EV_Cell_Gene_ID_moduleColor_sft_pwr_14_h0.2, file = "rna_WGCNA_EV_Cell_Gene_ID_moduleColor_sft_pwr_14_h0.2.csv")


collectGarbage()

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute"))
BiocManager::install("WGCNA")
install.packages("BiocManager")
library(WGCNA)
library(flashClust)
library(cluster)
library(DOSE)
library(gplots)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#4 Interfacing network analysis with other data such as functional annotation and gene ontology
#Our previous analysis has identified several modules (labeled brown, red, and salmon) that are highly associated with weight. 
#To facilitate a biological interpretation, we would like to know the gene ontologies of the genes in the modules, 
#whether they are significantly enriched in certain functional categories etc.

#4.a Output gene lists for use with online software and services

#One option is to simply export a list of gene identifiers that can be used as input for several popular gene ontology
#and functional enrichment analysis suites such as David or AmiGO. 
#For example, we write out the LocusLinkID (entrez) codes for the brown module into a file:

View(datExpr)
# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExpr)
#probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = geneInfo$ENTREZID;
# $ Choose interesting modules
intModules = c("blue", "brown", "yellow", "green", "turquoise")
for (module in intModules)
{
  # Select module probes
  modGenes = (dynamicColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
#probes = names(datExpr)
#probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
#allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
#intModules = c("brown", "red", "salmon")
#for (module in intModules)
#{
# Select module probes
#  modGenes = (moduleColors==module)
# Get their entrez ID codes
#  modLLIDs = allLLIDs[modGenes];
# Write them into a file
#  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
#  write.table(as.data.frame(modLLIDs), file = fileName,
#              row.names = FALSE, col.names = FALSE)
#}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)
View(allLLIDs)


GOenr = GOenrichmentAnalysis(mergedColors, allLLIDs, organism = "mouse", nBestP = 10);

tab = GOenr$bestPTerms[[4]]$enrichment

View(tab)
names(tab)

write.table(tab, file = "rna_WGCNA_EV_Cell_GOEnrichmentTable_sft_pwer_14_h0.2.csv", sep = ",", quote = TRUE, row.names = FALSE)



keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab
View(screenTab)
View(tab)
write.table(screenTab, file = "rna_WGCNA_EV_Cell_GOEnrichmentTable_shortened_sft_pwr_14_h0.2.csv", sep = ",", quote = TRUE, row.names = FALSE)

save.image("Transcriptomic_WGCNA_v4.RData")
getwd()

####Consensus Modules####
###To compare and contrast EV vs Cell groups, let's find consensus modules (retained modules in both data/groups)
#for EV and Cell groups, look at the genes and pathways,
#then, maybe look at modules are most different between those two groups. 

getwd()
#install stuff first
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute"))
BiocManager::install("WGCNA")
install.packages("BiocManager")
install.packages("Seurat")
library(WGCNA)
library(flashClust)
library(cluster)
library(DOSE)
library(gplots)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()
#practice first 
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-DataInput.pdf
#From horvath's tutorial:
#II. Consensus network analysis of liver expression data, female and male mice
#1. Data input and cleaning

#1.a Loading expression data
#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv");
# Read in the male liver data set
maleData = read.csv("LiverMale3600.csv");
# Take a quick look at what is in the data sets (caution, longish output):
dim(femData)
names(femData)
dim(maleData)
names(maleData)
View(femData)

# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(femData[-c(1:8)])));
names(multiExpr[[1]]$data) = femData$substanceBXH;
rownames(multiExpr[[1]]$data) = names(femData)[-c(1:8)];
multiExpr[[2]] = list(data = as.data.frame(t(maleData[-c(1:8)])));
names(multiExpr[[2]]$data) = maleData$substanceBXH;
rownames(multiExpr[[2]]$data) = names(maleData)[-c(1:8)];
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)
exprSize
#$nSets
#[1] 2

#$nGenes
#[1] 3600

#$nSamples
#[1] 135 124

#$structureOK
#[1] TRUE

#1.b Rudimentary data cleaning and outlier removal
#We first identify genes and samples with excessive numbers of missing samples. These can be identified using the
#function goodSamplesGenesMS.

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

#If the last statement returns TRUE, all genes and samples have passed the cuts. If it returns FALSE, the following
#code removes the offending samples and genes:
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

#We now cluster the samples on their Euclidean distance, separately in each set.
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
#The easiest way to see the two dendrograms at the same time is to plot both into a pdf file that can be viewed
#using standard pdf viewers.
pdf(file = "Plots_SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();

#By inspection, there seems to be one outlier in the female data set, and no obvious outliers in the male set. We
#now remove the female outlier using a semi-automatic code that only requires a choice of a height cut. We first
#re-plot the two sample trees with the cut lines included (see Fig. 1):

# Choose the "base" cut height for the female data set
baseHeight = 16
# Adjust the cut height for the male data set for the number of samples
cutHeights = c(16, 16*exprSize$nSamples[2]/exprSize$nSamples[1]);
# Re-plot the dendrograms including the cut lines
pdf(file = "Plots_SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
{
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
  abline(h=cutHeights[set], col = "red");
}
dev.off();

#Now comes the actual outlier removal:
for (set in 1:nSets)
{
  # Find clusters cut by the line
  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
  # Keep the largest one (labeled by the number 1)
  keep = (labels==1)
  multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
}

collectGarbage();
# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize

#$nSets
#[1] 2

#$nGenes
#[1] 3600

#$nSamples
#[1] 134 124

#$structureOK
#[1] TRUE

#we can see it removed that one female outlier sample

#1.c Loading clinical trait data
#We now read in the trait data and match the samples for which they were measured to the expression samples.
traitData = read.csv("ClinicalTraits.csv");
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];
# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)
allTraits$Mice
# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, allTraits$Mice);
  Traits[[set]] = list(data = allTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;
#We now have the expression data for both sets in the variable multiExpr, and the corresponding clinical traits in
#the variable Traits. The last step is to save the relevant data for use in the subsequent analysis:
save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,
     file = "Consensus-dataInput.RData");

#2.a One-step automatic network construction and module detection
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData");

# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets


#Among other variables we have loaded the variables multiExpr and Traits containing the expression and trait
#data, respectively. Further, expression data dimensions are stored in nGenes and nSamples.
#2 Network construction and module detection
#This step is the bedrock of all network analyses using the WGCNA methodology. We present three different ways
#of constructing a network and identifying modules:
#  a. Using a convenient 1-step function for network construction and detection of consensus modules, suitable for
#users wishing to arrive at the result with minimum effort;
#b. Step-by-step network construction and module detection for users who would like to experiment with customized/alternate methods;
#c. An automatic block-wise network construction and module detection method for users who wish to analyze
#data sets too large to be analyzed all in one.
#In this tutorial section, we illustrate the 1-step, automatic multiple set network construction and detection of
#consensus modules. We note that while the actual network construction and module detection is executed in a
#single function call, a preliminary step of choosing a suitable soft-thresholding power must be performed first.

#for Dae's data, we will use "b". "a" is fine too. "c" is to save memory space, but I won't be using that.

#2.a One-step network construction and module detection
#2.a.1 Choosing the soft-thresholding power: analysis of network topology
#Constructing a weighted gene network entails the choice of the soft thresholding power β to which co-expression
#similarity is raised to calculate adjacency [1]. The authors of [1] have proposed to choose the soft thresholding
#power based on the criterion of approximate scale-free topology. We refer the reader to that work for more details;
#here we illustrate the use of the function pickSoftThreshold that performs the analysis of network topology and
#aids the user in choosing a proper soft-thresholding power. The user chooses a set of candidate powers (the function
#provides suitable default values), and the function returns a set of network indices that should be inspected.

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)

colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}

#The result is shown in Plots. We choose the power 6 for both sets, 
#which is the first power that is over the 0.8 threshold.

#2.a.2 Calculation of network adjacencies
#Network construction starts by calculating the adjacencies in the individual sets, 
#using the soft thresholding power 6:
softPower = 6;
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower;

#2.a.3 Calculation of Topological Overlap
#We now turn the adjacencies into Topological Overlap Matrix (TOM) [2, 3]:

# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ]);

#2.a.4 Scaling of Topological Overlap Matrices to make them comparable across sets
#Topological Overlap Matrices of different data sets may have different statistical properties. For example, the
#TOM in the male data may be systematically lower than the TOM in female data. Since consensus is defined
#as the component-wise minimum of the two TOMs, a bias may result. Here we illustrate a simple scaling that
#mitigates the effect of different statistical properties to some degree. We scale the male TOM such that the 95th
#percentile equals the 95th percentile of the female TOM:

# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(081722)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets)
{
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
  # Scale the male TOM
  if (set>1)
  {
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
    TOM[set, ,] = TOM[set, ,]^scalePowers[set];
  }
}

#The array TOM now contains the scaled TOMs. To see what the scaling achieved, we form a quantile-quantile plot
#of the male and female topological overlaps before and after scaling:

# For plotting, also scale the sampled TOM entries
scaledTOMSamples = list();
for (set in 1:nSets)
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
# Open a suitably sized graphics window
sizeGrWindow(6,6)
pdf(file = "Plots_TOMScaling-QQPlot.pdf", wi = 6, he = 6);
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off();

#The result is shown in Fig. 2. 
#In this case the scaling changed the male TOM only very slightly, and brought it
#closer to the reference line shown in blue.

#2.a.5 Calculation of consensus Topological Overlap
#We now calculate the consensus Topological Overlap by taking the component-wise (“parallel”) minimum of the
#TOMs in individual sets    
consensusTOM = pmin(TOM[1, , ], TOM[2, , ]);
#Thus, the consensus topological overlap of two genes is only large if the corresponding entries in the two sets are
#also large.

#2.a.6 Clustering and module identification
#We use the consensus TOM as input to hierarchical clustering, and identify modules in the resulting dendrogram
#using the Dynamic Tree Cut algorithm [1]
# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average");
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)
#To see a quick summary of the module detection, we use table(unmergedLabels):
table(unmergedLabels)
#unmergedLabels
#0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19 
#351 352 300 285 281 264 259 244 239 184 127 127 105 103  90  82  64  62  44  37 
#The Dynamic Tree Cut returned 19 proper modules with sizes ranging from 352 to 37 genes. The label 0 is
#reserved for genes not assigned to any of the modules. The following code plots the consensus gene dendrogram
#together with the preliminary module colors:
sizeGrWindow(8,6)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#The resulting plot is shown in Plots.

#2.a.7 Merging of modules whose expression profiles are very similar
#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent
#to merge such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire
#modules, we calculate their eigengenes (MEs) and cluster them on their consensus correlation, that is the minimum
#correlation across the two sets:

## Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");
# Plot the result
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col = "red")
#The resulting tree is shown in Fig. 4. 
#Figure 4: Dendrogram of consensus module eigengenes obtained by clustering the eigengenes on their consensus
#correlation across male and female mice. The red line is the merging threshold; groups of eigengenes below the
#threshold represent modules whose expression profiles are too similar and should be merged.
#Two pairs of modules falls below the merging threshold. The merging can
#be performed automatically:
merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)  
#The variable merge contains various information; we will need the following:
# Numeric module labels
moduleLabels = merge$colors;
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs;
#Lastly, we plot the gene dendrogram again, this time with both the unmerged and the merged module colors:
sizeGrWindow(9,6)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#The plot is shown in Fig. 5. We now save the information necessary in the subsequent parts of the tutorial:
save(consMEs, moduleColors, moduleLabels, consTree, file = "Consensus-NetworkConstruction-man.RData")
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")



###let's try up to here with Dae's data####
library(WGCNA)
library(flashClust)
library(cluster)
library(DOSE)
library(gplots)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()
#practice first 
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-DataInput.pdf
#From horvath's tutorial:
#II. Consensus network analysis of liver expression data, female and male mice
#1. Data input and cleaning

#1.a Loading expression data
#Read in the EV data set
View(datExprEV)
# Read in the Cell data set
View(datExprCell)
# Take a quick look at what is in the data sets (caution, longish output):
dim(datExprEV)
names(datExprEV)
dim(datExprCell)
names(datExprCell)

View(femData)
# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("EV samples", "Cell samples")
shortLabels = c("EV", "Cell")
# Form multi-set expression data: tutorial's femData has columns starting from 9 contain actual expression data,
#with genes as rows and samples as columns.
#My samples have samples as rows and genes as columns.


multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame((datExprEV)));
names(multiExpr[[1]]$data) = colnames(datExprEV);
rownames(multiExpr[[1]]$data) = colnames(t(datExprEV));
multiExpr[[2]] = list(data = as.data.frame((datExprCell)));
names(multiExpr[[2]]$data) = colnames(datExprCell);
rownames(multiExpr[[2]]$data) = colnames(t(datExprCell));

#names(multiExpr[[1]]$data)
#rownames(multiExpr[[1]]$data)


# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)
exprSize
#$nSets
#[1] 2

#$nGenes
#[1] 19768

#$nSamples
#[1] 18 18

#$structureOK
#[1] TRUE

#1.b Rudimentary data cleaning and outlier removal
#We first identify genes and samples with excessive numbers of missing samples. These can be identified using the
#function goodSamplesGenesMS.

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

#If the last statement returns TRUE, all genes and samples have passed the cuts. If it returns FALSE, the following
#code removes the offending samples and genes:
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

#We now cluster the samples on their Euclidean distance, separately in each set.
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
#The easiest way to see the two dendrograms at the same time is to plot both into a pdf file that can be viewed
#using standard pdf viewers.
pdf(file = "Plots_SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();

#This plot is to find sample outlier, but since I only have 36 samples, let's not remove any. 


#1.c Loading clinical trait data
#We now read in the trait data and match the samples for which they were measured to the expression samples.
View(traitData)
View(allTraits)
View(datTraits)
#datTraits is dae's
#let's create another column with sample name
allTraits<-datTraits
allTraits$Sample<-rownames(allTraits)
allTraits<-allTraits[,c(5,1:4)]

View(Traits)
# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)
allTraits$Sample

# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, allTraits$Sample);
  Traits[[set]] = list(data = allTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;
#We now have the expression data for both sets in the variable multiExpr, and the corresponding clinical traits in
#the variable Traits. The last step is to save the relevant data for use in the subsequent analysis:
save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,
     file = "dae_Consensus-dataInput.RData");

getwd()
save.image("Transcriptomic_WCGNA_v5.RData")


#2.a One-step automatic network construction and module detection
# Load the data saved in the first part
lnames = load(file = "dae_Consensus-dataInput.RData");

# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets


#Among other variables we have loaded the variables multiExpr and Traits containing the expression and trait
#data, respectively. Further, expression data dimensions are stored in nGenes and nSamples.
#2 Network construction and module detection
#This step is the bedrock of all network analyses using the WGCNA methodology. We present three different ways
#of constructing a network and identifying modules:
#  a. Using a convenient 1-step function for network construction and detection of consensus modules, suitable for
#users wishing to arrive at the result with minimum effort;
#b. Step-by-step network construction and module detection for users who would like to experiment with customized/alternate methods;
#c. An automatic block-wise network construction and module detection method for users who wish to analyze
#data sets too large to be analyzed all in one.
#In this tutorial section, we illustrate the 1-step, automatic multiple set network construction and detection of
#consensus modules. We note that while the actual network construction and module detection is executed in a
#single function call, a preliminary step of choosing a suitable soft-thresholding power must be performed first.

#for Dae's data, we will use "b". "a" is fine too. "c" is to save memory space, but I won't be using that.

#2.a One-step network construction and module detection
#2.a.1 Choosing the soft-thresholding power: analysis of network topology
#Constructing a weighted gene network entails the choice of the soft thresholding power β to which co-expression
#similarity is raised to calculate adjacency [1]. The authors of [1] have proposed to choose the soft thresholding
#power based on the criterion of approximate scale-free topology. We refer the reader to that work for more details;
#here we illustrate the use of the function pickSoftThreshold that performs the analysis of network topology and
#aids the user in choosing a proper soft-thresholding power. The user chooses a set of candidate powers (the function
#provides suitable default values), and the function returns a set of network indices that should be inspected.

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)

colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}

#The result is shown in Plots. We choose the power 4 for both sets, 
#which is the first power that is over the 0.8 threshold.


#2.a.2 Calculation of network adjacencies
#Network construction starts by calculating the adjacencies in the individual sets, 
#using the soft thresholding power 4:
softPower = 4;
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower;

#2.a.3 Calculation of Topological Overlap
#We now turn the adjacencies into Topological Overlap Matrix (TOM) [2, 3]:

# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ]);

#2.a.4 Scaling of Topological Overlap Matrices to make them comparable across sets
#Topological Overlap Matrices of different data sets may have different statistical properties. For example, the
#TOM in the male data may be systematically lower than the TOM in female data. Since consensus is defined
#as the component-wise minimum of the two TOMs, a bias may result. Here we illustrate a simple scaling that
#mitigates the effect of different statistical properties to some degree. We scale the EV TOM such that the 95th
#percentile equals the 95th percentile of the Cell TOM:

# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(081922)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets)
{
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
  # Scale the male TOM
  if (set>1)
  {
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
    TOM[set, ,] = TOM[set, ,]^scalePowers[set];
  }
}

#The array TOM now contains the scaled TOMs. To see what the scaling achieved, we form a quantile-quantile plot
#of the Ev and Cell topological overlaps before and after scaling:

# For plotting, also scale the sampled TOM entries
scaledTOMSamples = list();
for (set in 1:nSets)
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
# Open a suitably sized graphics window
sizeGrWindow(6,6)
pdf(file = "Plots_TOMScaling-QQPlot.pdf", wi = 6, he = 6);
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off();

#The result is shown in Fig. 2. 
#In this case the scaling changed the Cell TOM very slightly, and brought it
#closer to the reference line shown in blue.

#2.a.5 Calculation of consensus Topological Overlap
#We now calculate the consensus Topological Overlap by taking the component-wise (“parallel”) minimum of the
#TOMs in individual sets    
consensusTOM = pmin(TOM[1, , ], TOM[2, , ]);
#Thus, the consensus topological overlap of two genes is only large if the corresponding entries in the two sets are
#also large.


#2.a.6 Clustering and module identification
#We use the consensus TOM as input to hierarchical clustering, and identify modules in the resulting dendrogram
#using the Dynamic Tree Cut algorithm [1]
# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average");
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)
#To see a quick summary of the module detection, we use table(unmergedLabels):
table(unmergedLabels)
#unmergedLabels
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21 
#3063 2202 1941 1241  871  496  340  339  320  311  294  260  250  235  233  218  216  205  204  203  201 
#22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42 
#197  195  189  185  182  178  175  171  167  165  156  156  147  147  146  146  145  139  138  137  135 
#43   44   45   46   47   48   49   50   51   52   53   54   55   56   57   58   59   60   61   62   63 
#128  126  123  120  115  114  108  108   98   91   85   83   82   80   79   78   78   73   72   67   66 
#64   65   66   67   68   69   70   71   72   73   74   75   76   77 
#63   61   61   59   55   55   55   54   53   53   52   46   46   42 

#tutorial one gave 19
#0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19 
#351 352 300 285 281 264 259 244 239 184 127 127 105 103  90  82  64  62  44  37 
#The Dynamic Tree Cut returned 77 proper modules with sizes ranging from 3063 to 42 genes. 
#The label 0 is reserved for genes not assigned to any of the modules. 
#The following code plots the consensus gene dendrogram
#together with the preliminary module colors:
sizeGrWindow(8,6)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#The resulting plot is shown in Plots.

#2.a.7 Merging of modules whose expression profiles are very similar
#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent
#to merge such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire
#modules, we calculate their eigengenes (MEs) and cluster them on their consensus correlation, that is the minimum
#correlation across the two sets:

## Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average");
# Plot the result
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col = "red")
#The resulting tree is shown in Fig. 4. 
#Figure 4: Dendrogram of consensus module eigengenes obtained by clustering the eigengenes on their consensus
#correlation across EV and Cell mice. The red line is the merging threshold; groups of eigengenes below the
#threshold represent modules whose expression profiles are too similar and should be merged.
#Two pairs of modules falls below the merging threshold. The merging can
#be performed automatically:
merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)  
#The variable merge contains various information; we will need the following:
# Numeric module labels
moduleLabels = merge$colors;
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs;
#Lastly, we plot the gene dendrogram again, this time with both the unmerged and the merged module colors:
sizeGrWindow(9,6)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#The plot is shown in Fig. 5. We now save the information necessary in the subsequent parts of the tutorial:
save(consMEs, moduleColors, moduleLabels, consTree, file = "dae_Consensus-NetworkConstruction-man.RData")
save.image("Transcriptomic_WCGNA_v4.RData")




####Consensus Modules continued####



library(WGCNA)
library(flashClust)
library(cluster)
library(DOSE)
library(gplots)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-RelateToFemMods.pdf  
#II. Consensus network analysis of expression data, ev and cell data sets
#3. Relating consensus modules to female set-specific modules  
#0 Preliminaries: setting up the R session
#Here we assume that a new R session has just been started. We load the WGCNA package, set up basic parameters
#and load data saved in the parts 1 and 2 of the tutorial.
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.

#setwd();
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the data saved in the first part
lnames = load(file = "dae_Consensus-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load the results of network analysis, tutorial part 2.a
lnames = load(file = "dae_Consensus-NetworkConstruction-man.RData");
lnames

#We have loaded the variables multiExpr and Traits containing the expression and trait data, respectively. Further,
#expression data dimensions are stored in nGenes and nSamples.

#now for ev-specific vs consensus####
#3 Relating consensus modules to ev set-specific modules
#In this section we assume the reader has worked through the one-step network construction and module detection
#in the ev data (Section 2.a of the female mouse liver analysis tutorial). We first load the ev
#data and rename them so that they do not conflict with the consensus data:

#save(MEs, dynamicModsEV, dynamicColorsEV, geneTreeEV,
#     file = "EV-02-networkConstruction-auto.RData")
#save(mergedMEsEV, moduleLabelsEV, moduleColorsEV, geneTreeEV,
#file = "EV-02-networkConstruction-man.RData")

#lnames = load("EV-02-networkConstruction-auto.RData")
#lnames

lnames = load("EV-02-networkConstruction-man.RData")
lnames

# Rename variables to avoid conflicts
EVLabels = moduleLabelsEV;
EVColors = moduleColorsEV;
EVTree = geneTreeEV;
EVMEs = orderMEs(mergedMEsEV, greyName = "ME0");

#Next we load the results of the consensus module identification:
lnames = load("dae_Consensus-NetworkConstruction-man.RData")
lnames

#The consensus network analysis results are represented by the variables consMEs, moduleLabels, moduleColors, and
#consTree. We are now ready to relate the EV modules to the consensus modules. We calculate the overlaps
#of each pair of EV-consensus modules, and use the Fisher’s exact test (also known as hypergeometric test) to
#assign a p-value to each of the pairwise overlaps.
#Isolate the module labels in the order they appear in ordered module eigengenes
EVModuleLabels = substring(names(EVMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
EVModules = labels2colors(as.numeric(EVModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))
# Numbers of female and consensus modules
nEVMods = length(EVModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nEVMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nEVMods, ncol = nConsMods);
# Execute all pairwise comparisons
for (EVmod in 1:nEVMods)
  for (cmod in 1:nConsMods)
  {
    EVMembers = (EVColors == EVModules[EVmod]);
    consMembers = (moduleColors == consModules[cmod]);
    pTable[EVmod, cmod] = -log10(fisher.test(EVMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[EVmod, cmod] = sum(EVColors == EVModules[EVmod] & moduleColors ==
                                  consModules[cmod])
  }
#To display the p-value and count tables in an informative way, we create a color-coded table of the intersection
#counts. The colors will indicate the p-value significance:

# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
EVModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# Actual plotting
sizeGrWindow(10,7 );
pdf(file = "Plots_ConsensusVsEVModules.pdf", wi = 20, he = 14);
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", EVModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("EV ", EVModules, ": ", EVModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of EV sample-specific and EV-Cell consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
dev.off();



#The resulting color-coded table is shown in Fig. 1. The table indicates that most female set-specific modules have
#a consensus counterpart. This indirectly shows that the module structure in the male liver expression data is
#very similar to the female data. Interestingly, there are two female modules (labeled by grey60 and lightgreen
#colors) that have no consensus counterpart; almost all genes in the two female modules are labeled “grey”, that
#is unassigned, in the consensus network.


#now for cell-specific vs consensus####

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the data saved in the first part
lnames = load(file = "dae_Consensus-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load the results of network analysis, tutorial part 2.a
lnames = load(file = "dae_Consensus-NetworkConstruction-man.RData");
lnames

#We have loaded the variables multiExpr and Traits containing the expression and trait data, respectively. Further,
#expression data dimensions are stored in nGenes and nSamples.
#3 Relating consensus modules to female set-specific modules
#In this section we assume the reader has worked through the one-step network construction and module detection
#in the female mouse liver data (Section 2.a of the female mouse liver analysis tutorial). We first load the female
#data and rename them so that they do not conflict with the consensus data:

lnames = load("Cell-02-networkConstruction-man.RData")
lnames
#"mergedMEsCell"    "moduleLabelsCell" "moduleColorsCell" "geneTreeCell"    

# Rename variables to avoid conflicts
CellLabels = moduleLabelsCell;
CellColors = moduleColorsCell;
CellTree = geneTreeCell;
CellMEs = orderMEs(mergedMEsCell, greyName = "ME0");

#Next we load the results of the consensus module identification:
lnames = load("dae_Consensus-NetworkConstruction-man.RData")
lnames
#"consMEs"      "moduleColors" "moduleLabels" "consTree"    


#The consensus network analysis results are represented by the variables consMEs, moduleLabels, moduleColors, and
#consTree. We are now ready to relate the Cell modules to the consensus modules. We calculate the overlaps
#of each pair of Cell-consensus modules, and use the Fisher’s exact test (also known as hypergeometric test) to
#assign a p-value to each of the pairwise overlaps.
#Isolate the module labels in the order they appear in ordered module eigengenes
CellModuleLabels = substring(names(CellMEs), 3)
consModuleLabels = substring(names(consMEs[[2]]$data), 3)
#consMEs[[1]] are 18 EV samples and consMEs[[2]] are for 18 Cell samples,
#but when you are doing names(consMEs[[1]]$data), then names(consMEs[[1]]$data) or names(consMEs[[2]$data) are same

# Convert the numeric module labels to color labels
CellModules = labels2colors(as.numeric(CellModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))
# Numbers of female and consensus modules
nCellMods = length(CellModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nCellMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nCellMods, ncol = nConsMods);
# Execute all pairwise comparisons
for (Cellmod in 1:nCellMods)
  for (cmod in 1:nConsMods)
  {
    CellMembers = (CellColors == CellModules[Cellmod]);
    consMembers = (moduleColors == consModules[cmod]);
    pTable[Cellmod, cmod] = -log10(fisher.test(CellMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[Cellmod, cmod] = sum(CellColors == CellModules[Cellmod] & moduleColors ==
                                    consModules[cmod])
  }
#To display the p-value and count tables in an informative way, we create a color-coded table of the intersection
#counts. The colors will indicate the p-value significance:

# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
CellModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# Actual plotting
sizeGrWindow(10,7 );
pdf(file = "Plots_ConsensusVsCellModules.pdf", wi = 20, he = 14);
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", CellModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("Cell ", CellModules, ": ", CellModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of Cell sample-specific and EV-Cell consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
dev.off();



#The resulting color-coded table is shown in Fig. 1. The table indicates that most female set-specific modules have
#a consensus counterpart. This indirectly shows that the module structure in the male liver expression data is
#very similar to the female data. Interestingly, there are two female modules (labeled by grey60 and lightgreen
#colors) that have no consensus counterpart; almost all genes in the two female modules are labeled “grey”, that
#is unassigned, in the consensus network.

#let's do correspondence of ev-specific modules vs cell-specific modules
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nEVMods, ncol = nCellMods);
CountTbl = matrix(0, nrow = nEVMods, ncol = nCellMods);
# Execute all pairwise comparisons
for (EVmod in 1:nEVMods)
  for (Cellmod in 1:nCellMods)
  {
    EVMembers = (EVColors == EVModules[EVmod]);
    CellMembers = (CellColors == CellModules[Cellmod]);
    pTable[EVmod, Cellmod] = -log10(fisher.test(EVMembers, CellMembers, alternative = "greater")$p.value);
    CountTbl[EVmod, Cellmod] = sum(EVColors == EVModules[EVmod] & CellColors ==
                                     CellModules[Cellmod])
  }
#To display the p-value and count tables in an informative way, we create a color-coded table of the intersection
#counts. The colors will indicate the p-value significance:

# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;

sizeGrWindow(10,7 );
pdf(file = "Plots_EVModulesVsCellModules.pdf", wi = 20, he = 14);
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", CellModules),
               yLabels = paste(" ", EVModules),
               colorLabels = TRUE,
               xSymbols = paste("Cell ", CellModules, ": ", CellModTotals, sep=""),
               ySymbols = paste("EV ", EVModules, ": ", EVModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of EV sample-specific and Cell sample-specific modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
dev.off();

getwd()

save.image("Transcriptomic_WCGNA_v5.RData")


###  4. Relating consensus modules to external microarray sample information and exporting network analysis results####
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-RelateModsToTraits.pdf

#In this section we illustrate the use of module eigengenes to relate consensus modules to external microarray
#sample information such as classical clinical traits. In this analysis we have available several clinical traits. We
#relate the traits to consensus module eigengenes in each of the two sets. We remind the reader that while the
#consensus modules is a single module assignment for all genes, the module eigengenes represent the modules in
#each of the two sets. In other words, we have a single module assignment for each gene, but we have two sets
#of consensus module eigengenes, because a given module (set of genes) has a particular expression profile in the
#female mice, and a different expression profile in the male mice. Similarly, we have the trait data separately for
#the female and for the male mice.

#Here we assume that a new R session has just been started. We load the WGCNA package, set up basic parameters
#and load data saved in the parts 1 and 2 of the tutorial.
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the data saved in the first part
lnames = load(file = "dae_Consensus-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
#[1] "multiExpr"   "Traits"      "nGenes"      "nSamples"    "setLabels"   "shortLabels" "exprSize" 
# Also load results of network analysis
lnames = load(file = "dae_Consensus-NetworkConstruction-man.RData");
lnames
#[1] "consMEs"      "moduleColors" "moduleLabels" "consTree"    

exprSize = checkSets(multiExpr);
nSets = exprSize$nSets;

#We have loaded the variables multiExpr and Traits containing the expression and trait data, respectively. Further, expression data dimensions are stored in nGenes and nSamples. The consensus network analysis results are
#represented by the variables consMEs, moduleLabels, moduleColors, and consTree.

# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in 1:nSets)
{
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}

#We now display the module-trait relationships using a color-coded table. We print the correlations and the
#corresponding p-values, and color-code the entris by the p-value significance.

# Convert numerical labels to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
# Open a suitably sized window (the user should change the window size if necessary)
sizeGrWindow(10,7)
pdf(file = "Plots_Consensus_Module-trait_relationship_in_EV_Samples.pdf", wi = 20, he = 14);
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus Module-trait relationships in", setLabels[set]))
dev.off();
# Plot the module-trait relationship table for set number 2

pdf(file = "Plots_Consensus_Module-trait_relationship_in_Cell_Samples.pdf", wi = 20, he = 14);
set = 2
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus Module-trait relationships in", setLabels[set]))
dev.off();

#The two tables are shown in Figs. 1 and 2. 
#Figure 1: Relationships of consensus module eigengenes and clinical traits in the female data. Each row in the
#table corresponds to a consensus module, and each column to a trait. Numbers in the table report the correlations
#of the corresponding module eigengenes and traits, with the p-values printed below the correlations in parentheses.
#The table is color coded by correlation according to the color legend.

#The two module-trait relationship tables look similar but they are
#not the same. For example, they both identify the turquoise, purple and green modules as highly related to
#weight, although the actual correlations and p-values differ slightly. There are several ways of forming a measure
#of module-trait relationships that summarizes the two sets into one measure. We will form a very conservative
#one: for each module-trait pair we take the correlation that has the lower absolute value in the two sets if the
#two correlations have the same sign, and zero relationship if the two correlations have opposite signs:

#(side note, because Traits has EV vs cell, pregnancy, sex, and ethanol, 
#there is absolutely no consensus correlation or p-value for EV vs cell Trait, 
#because the correlations in the EV and cell data sets will all have opposite signs for trait Location (EV vs cell)
#and no consensus can be formed.
#therefore, because we will end up getting an error
#Error in consensusCor[negative] = pmax(moduleTraitCor[[1]][negative],  : NAs are not allowed in subscripted assignments
#we have to exclude the trait Location (EV vs cell) to make sure the codes work.)
#rows are MEs values and columns are traits with the first column being Location
#so we will exclude the first column using: moduleTraitCor[[1]][,c(2:4)] instead of moduleTraitCor[[1]]
#pretty much, i am adding [,c(2:4)] to bunch of stuff
#pretty easy, but had to figure out names(Traits[[set]]$data[,c(2:4)]) to solve another error

# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]][,c(2:4)]), ncol(moduleTraitCor[[1]][,c(2:4)]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]][,c(2:4)]), ncol(moduleTraitCor[[1]][,c(2:4)]));
# Find consensus negative correlations
negative = moduleTraitCor[[1]][,c(2:4)] < 0 & moduleTraitCor[[2]][,c(2:4)] < 0;
consensusCor[negative] = pmax(moduleTraitCor[[1]][,c(2:4)][negative], moduleTraitCor[[2]][,c(2:4)][negative]);
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][,c(2:4)][negative], moduleTraitPvalue[[2]][,c(2:4)][negative]);
# Find consensus positive correlations
positive = moduleTraitCor[[1]][,c(2:4)] > 0 & moduleTraitCor[[2]][,c(2:4)] > 0;
consensusCor[positive] = pmin(moduleTraitCor[[1]][,c(2:4)][positive], moduleTraitCor[[2]][,c(2:4)][positive]);
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][,c(2:4)][positive], moduleTraitPvalue[[2]][,c(2:4)][positive]);

#We display the consensus module–trait relationships again using a color-coded table:

sizeGrWindow(10,7)
pdf(file = "Plots_ModuleTraitRelationships-consensus.pdf", wi = 20, he = 14);
textMatrix = paste(signif(consensusCor, 2), "\n(",
                   signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]][,c(2:4)])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]$data[,c(2:4)]),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))
dev.off();

#The table is shown in Fig. 3. 
#Figure 3: Consensus relationships of consensus module eigengenes and clinical traits across the female and male
#data. Each row in the table corresponds to a consensus module, and each column to a trait. Numbers in the table
#report the consensus correlations of the corresponding module eigengenes and traits, with the p-values printed
#below the correlations in parentheses. The table is color coded by correlation according to the color legend.
#Missing (NA) entries indicate that the correlations in the male and female data sets have opposite signs and no
#consensus can be formed.

#The advantage of the consensus relationship table is that it isolates the module-trait
#relationships that are present in both sets, and hence may be in a sense considered validated. For example, we
#confirm that the turquoise, purple, and green modules are highly related to weight in both sets; the brown module
#is highly related to insulin levels etc. One could now adapt the gene selection technique illustrated in the female
#expression analysis tutorial to look for particular genes that are highly related to traits as well as being highly
#connected in a module related to the trait; we leave this analysis as an exercise for the reader.

#4.a Exporting results of the network analysis
#We now put together a data frame that summarizes the results of network analysis, namely the gene significances (GS) 
#and module memberships (also known as kME) of all probes. We start by loading the gene annotation table.

View(annot)
View(datExpr)
View(tdatExpr)
View(tdatExpr_Entrez)

#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
probes = names(multiExpr[[1]]$data)
probes2annot = match(probes, tdatExpr$Gene_ID)

#We next (re-)calculate the module eigengenes in the ”alphabetic” order and calculate the gene significances and
#module memberships in each data set.
consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list();
kME = list();
for (set in 1:nSets)
{
  GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data);
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}

#We perform a very simple ”meta-analysis” by combining the Z scores of correlations from each set to form a
#meta-Z score and the corresponding p-value.
GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z)/sqrt(2);
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2);
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);

#Next we form matrices holding the GS and kME. We use a simple re-shaping trick to put the values and the
#associated p-values and meta-analysis results next to one another.
GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP);
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 6*nTraits)
rownames(GSmat) = probes;
colnames(GSmat) = spaste(
  c("GS.setEV.", "GS.setCell.", "p.GS.setEV.", "p.GS.setCell.", "Z.GS.meta.", "p.GS.meta"),
  rep(traitNames, rep(6, nTraits)))
# Same code for kME:
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
  c("kME.setEV.", "kME.setCell.", "p.kME.setEV.", "p.kME.setCell.", "Z.kME.meta.", "p.kME.meta"),
  rep(MEnames, rep(6, nMEs)))

#Finally we put together the full information data frame and write it into a plain text CSV file that can be read by
#standard spreadsheet programs. Note that the probes are not sorted in any particular way; many sort orders are
#possible and we leave it to the reader to either modify the code or to perform the sort in a spreadsheet software.
info = data.frame(Probe = probes, GeneSymbol = tdatExpr$Gene_ID[probes2annot],
                  EntrezID = tdatExpr_Entrez$ENTREZID[probes2annot],
                  ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels),
                  GSmat,
                  kMEmat);
write.csv(info, file = "rna_WGCNA_EV_Cell_consensusAnalysis-CombinedNetworkResults.csv",
          row.names = FALSE, quote = FALSE);

#5. Comparing eigengene networks in male and female mice
#In this section we compare the consensus eigengene networks in the female and male data sets (such comparison
#is often called differential analysis). Consensus eigengene networks capture the relationships among consensus
#modules; the relationships are quantified by eigengene correlations. We first extend the consensus eigengenes by
#adding the body weight clinical trait as an additional “eigengene”:

# Create a variable EV_Cell that will hold just the EV_Cell of mice in both sets
#DOESN'T work, because that sample trait EV_Cell will create NaN dissimilarity value, 
#since none of the EV samples will have Cell location and vice versa for Cell samples.
EV_Cell = vector(mode = "list", length = nSets);
for (set in 1:nSets)
{
  EV_Cell[[set]] = list(data = as.data.frame(Traits[[set]]$data$EV_Cell));
  names(EV_Cell[[set]]$data) = "EV_Cell"
}
# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors);
# We add the EV_Cell trait to the eigengenes and order them by consesus hierarchical clustering:
MET = consensusOrderMEs(addTraitToMEs(consMEsC, EV_Cell));


#So, Create a variable Sex that will hold just the Sex of mice in both sets
Sex = vector(mode = "list", length = nSets);
for (set in 1:nSets)
{
  Sex[[set]] = list(data = as.data.frame(Traits[[set]]$data$Sex));
  names(Sex[[set]]$data) = "Sex"
}
# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors);
# We add the weight trait to the eigengenes and order them by consesus hierarchical clustering:
MET = consensusOrderMEs(addTraitToMEs(consMEsC, Sex));

#We now call the function plotEigengeneNetworks that performs the differential analysis:
sizeGrWindow(8,10);
pdf(file = "Plots_EV_Cell_Consensus_EigengeneNetworks.pdf", wi = 20, he = 14);
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off();  

#Figure 1: Summary plot of consensus eigengene networks and their differential analysis. The top two panels
#show the dendrograms (clustering trees) of the consensus module eigengenes in the two sets indicated in the
#titles. Below, the eigengene networks in the two sets are shown as heatmaps labeled Female liver and Male liver.
#In the heatmaps, red denotes high adjacency (positive correlation) and blue denotes low adjacency (negative
#correlation). The Preservation heatmap shows the preservation network, defined as one minus the absolute
#difference of the eigengene networks in the two data sets. The barplot shows the mean preservation of adjacency
#for each of the eigengenes to all other eigengenes; in other words, the barplot depicts the column means of the
#preservation heatmap.

#The plot is shown in Fig. 1. The plot shows the eigengene dendrograms and eigengene network heatmaps,
#as well as two plots of network preservation between the two sets, namely a heatmap plot of the preservation
#adjacency [1] and a barplot of mean preservation of relationships for each eigengene. The overall preservation
#of the two eigengene networks is 0.71. A visual inspection of the EV and cell network
#heatmap plots indicates that the inter-module relationships in the two data sets are similar some modules and not similar in other modules. 
#The red blocks along the diagonal in the network heatmaps indicate meta-modules, groups of correlated eigengenes,
#and these are largely preserved between the two sets. The preservation heatmap and barplots indicate that some
#relationships are very highly preserved. The message here is that inter-module relationships are not strongly preserved
#but moderately preserved across EV and Cell data sets, meaning they may not be similar data sets,
#and encode biologically meaningful information.


#############################



###modules and traits####
getwd()
setwd("/scratch/user/daehyukchung")
#let's do a multidimensional scaling plot
#Multidimensional scaling (MDS) is a multivariate data analysis approach 
#that is used to visualize the similarity/dissimilarity between samples by plotting points in two dimensional plots.
cmd1EV=cmdscale(as.dist(dissTOMEV),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd1EV, col=as.character(colorh1_EV), main="EV MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")


######


#dynamicColors=staticColors
#below will give you module tables with gene names in each module
module_colorsEV= setdiff(unique(dynamicColorsEV), "grey")

for (color in module_colorsEV){
  moduleEV=SubGeneNamesEV[which(dynamicColorsEV==color)]
  write.table(moduleEV, paste("moduleEV_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}
#Dae, run colorh1_EV= dynamicColorsEV
#colorh1_EV= dynamicColorsEV
#Dae, to get hub genes, 
hubsEV    = chooseTopHubInEachModule(datExprEV, colorh1_EV)
hubsEV

#hubsEV
#black            blue           brown            cyan       darkgreen        darkgrey      darkorange 
#"Ttc29"         "Prmt1"        "Hs6st3"    "D8Ertd738e"          "Dpyd"        "Fam20a"          "Insc" 
#darkred   darkturquoise           green     greenyellow          grey60       lightcyan      lightgreen 
#"Bik"         "Prkcq"        "Anapc4"       "Nectin4"         "Atp5e" "1810007D17Rik"        "Lrrtm4" 
#lightyellow         magenta    midnightblue          orange            pink          purple             red 
#"Lhfpl3"         "Rtkn2"      "Serpini1"        "Gm1043"        "Ccser1"          "Cps1"       "Ccdc180" 
#royalblue          salmon         skyblue             tan       turquoise           white          yellow 
#"Mroh2b"       "Gm34256"       "Gm20471"        "Ifnar2"          "Scd2"         "Foxf2"          "Sfpq" 


####modules to external traits relationship###
#here, we want to see which modules relate closely to which traits
#for my samples, traits will be pregnancy/set, sex, and ethanol.
#it can also be ev or cell, if i have all the samples together.
#Module-Trait relationships. Color scale (red-green) represents the strength of the correlation between the module and the trait. 
#Each box gives a correlation value (R^2) followed by p-value (in parenthesis). 
#the color scale is the strength of the correlation between the module and the trait, 
#where for example, MEblack and Female being very red means highly correlated.

#let's load the excel csv file i have created for traits. 
#i am adding "row.names=1" function, which means make the first column (sample names) into row names 
#datTraits = read.csv("Traits_wgcna1.csv",row.names = 1)
#write.csv(datTraits,file="Traitswgcna.csv")
#form a data frame analogous to expression data that will hold the clinical traits.
#let's match the rownames of datExpr and datTraits, after you check that the samples are lined up correctly by rows
rownames(datTraits)
#let's get row 1:18 that only have EV samples, and column 2:4 to exclude EV_Cell trait
datTraitsEV<-datTraits[1:18,2:4]
rownames(datTraitsEV)
rownames(datExprEV)


rownames(datExprEV)<-rownames(datTraitsEV)

table(rownames(datTraitsEV)==rownames(datExprEV)) 
#should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraitsEV)

##
getwd()

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


#3 Relating modules to external clinical traits
#3.a Quantifying module-trait associations
#In this analysis we would like to identify modules that are significantly associated with the measured clinical traits.
#Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external
#traits and look for the most significant associations:

#instead of moduleColors, i have dynamicColorsEV from:
#dynamicColorsEV = labels2colors(dynamicModsEV)
View(dynamicColorsEV)
View(colorh1)

# Define numbers of genes and samples
nGenesEV = ncol(datExprEV);
nSamplesEV = nrow(datExprEV);
# Recalculate MEs with color labels
#moduleEigengenes() simply calculate the 1st Principal Component (PC) i.e., module eigengene (ME), of each module.
MEs0 = moduleEigengenes(datExprEV, dynamicColorsEV)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraitsEV, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamplesEV);

#let's use a graphical representation to read the table, a color-coded table,
#where we color code each association of module eigengenes to the sample trait by the correlation value:
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraitsEV),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("EV Module-Trait Relationships"))

#View(datTraits)
#colnames(datTraits)[1]<-"EV_Cell"
#colnames(datTraits)[2]<-"Pregnancy"
#colnames(datTraits)[4]<-"Ethanol"
#colnames(datTraits)[1]<-"EV_Cell"
#looking at module-trait relationships, there is a significant strong positive correlation of MEbrown, blue, turquoise modules to Cell trait,
#while rest are for EV trait. 


#Do barplot of trait-based mean gene significance across modules to find the module that has the highest gene significance to the trait. 
#III. Using simulated data to evaluate different module detection methods
#and gene screening approaches
#6. Relating modules and module eigengenes to external data
#6.a Representing modules by eigengenes and relating eigengenes to one another
#To get a sense of how related the modules are one can summarize each module by its eigengene (first principal component).


datMEEV=moduleEigengenes(datExprEV,colorh1_EV)$eigengenes
signif(cor(datMEEV, use="p"), 2)
#The result is



#We define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation
#between the module eigengenes, and use it to cluster the eigengene:
dissimMEEV=(1-t(cor(datMEEV, method="p")))/2
hclustdatMEEV=hclust(as.dist(dissimMEEV), method="average" )
# Plot the eigengene dendrogram
sizeGrWindow(8,9)
par(mfrow=c(1,1))
plot(hclustdatMEEV, main="EV Clustering tree based of the module eigengenes (first PC)")

hubsEV
hubsCell

#6.a.1 Pairwise scatter plots of the samples (arrays) along the module eigengenes
#We create a pairwise scatter plots of the samples (arrays) along the module eigengenes:
#define the microarray sample trait
y=datTraitsEV$Sex


sizeGrWindow(8,9)
plotMEpairs(datMEEV,y=y)
#Error in plot.new() : figure margins too large
#i got the error above, and all i had to do to fix it was to simply
#by making the sidebar bigger by clicking and dragging on its edge from right to left.
#dumb. 

#6.b Relating observed module eigengenes to simulated true module eigengenes
#We now relate the observed module eigengenes to the simulated true module eigengenes:
#attach(ModuleEigengeneNetwork1)
#ModuleEigengeneNetwork1=data.frame(y,MEturquoise,MEblue,MEbrown,MEyellow)

#6.e.2 Measure of module significance as average gene significance
#One can also define a measure of module significance as the average gene significance of all genes in the module. We
#use the absolute value for defining a correlation based gene significance measure.
GS1=as.numeric(cor(y,datExprEV, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, colorh1_EV, mean, na.rm=T)
ModuleSignificance
#here is the result:
#        black          blue         brown          cyan     darkgreen      darkgrey    darkorange       darkred darkturquoise 
#0.33197989    0.24618032    0.28281850    0.36751785    0.51035546    0.38973907    0.32463224    0.19662749    0.40150640 
#green   greenyellow          grey        grey60     lightcyan    lightgreen   lightyellow       magenta  midnightblue 
#0.28224869    0.20146168    0.52687155    0.35375722    0.15843617    0.45776639    0.34507046    0.28217802    0.25858798 
#orange          pink        purple           red     royalblue        salmon       skyblue           tan     turquoise 
#0.09472901    0.56497289    0.25756195    0.27099666    0.37365886    0.28062924    0.47649935    0.44819419    0.52859564 
#white        yellow 
#0.11382020    0.18181652 


#To plot module significance, one can use
sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,colorh1_EV)

#The advantage of this approach is that it can be used for any gene significance measure. 
#A gene significance measure could be defined without reference to a sample trait. 
#For example, it could indicate pathway membership (1 or 0) or gene essentiality (1 or 0), etc.
#The next logical step in an analysis of empirical data would be to carry out a functional enrichment analysis of
#each module, for example using the software EASE (David) at http://david.abcc.ncifcrf.gov/summary.jsp.
#We refer the reader to Tutorial I for an example of a functional enrichment analysis. Before we end, we again save
#calculated results for use in subsequent sections.
collectGarbage()
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")
getwd()

#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf

#3.b Gene relationship to trait and important modules: Gene Significance and Module Membership
#We quantify associations of individual genes with our trait of interest (Sex) by defining Gene Significance GS as
#(the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative
#measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 
#This allows us to quantify the similarity of all genes on the array to every module.
# Define variable trait containing the Sex column of datTrait
Sex = as.data.frame(datTraitsEV$Sex);
names(Sex) = "Sex"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExprEV, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExprEV, Sex, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Sex), sep="");
names(GSPvalue) = paste("p.GS.", names(Sex), sep="");


#3.c Intramodular analysis: identifying genes with high GS and MM
#Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module
#membership in interesting modules. As an example, we look at the lightcyan, magenta, purple, turquoise, blue module that has the highest association
#with Sex We plot a scatterplot of Gene Significance vs. Module Membership in the blue module:
#module lightcyon, magenta, purple highest gene significance for female, turquoise, blue for male.
# names (colors) of the modules
library(ggplot2)

module = "lightcyan"
column = match(module, modNames);
moduleGenes = dynamicColorsEV==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "cyan")

module = "magenta"
column = match(module, modNames);
moduleGenes = dynamicColorsEV==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "purple"
column = match(module, modNames);
moduleGenes = dynamicColorsEV==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module)


module = "turquoise"
column = match(module, modNames);
moduleGenes = dynamicColorsEV==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "blue"
column = match(module, modNames);
moduleGenes = dynamicColorsEV==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



#I. Network analysis of liver expression data in female mice
#5. Network visualization using WGCNA functions
#5.b Visualizing the network of eigengenes
#It is often interesting to study the relationships among the found modules. 
#One can use the eigengenes as representative profiles and quantify module similarity by eigengene correlation. 
#The package contains a convenient function plotEigengeneNetworks that generates a summary plot of the eigengene network. 
#It is usually informative to add a clinical trait (or multiple traits) to the eigengenes 
#to see how the traits fit into the eigengene network:
#Previously: 
#MEs0 = moduleEigengenes(datExprCell, dynamicColorsCell)$eigengenes
#MEs = orderMEs(MEs0)
## Isolate weight from the clinical traits
#Sex = as.data.frame(datTraits$Sex);
#names(Sex) = "Sex"
# Add the Sex to existing module eigengenes
MET = orderMEs(cbind(MEs, Sex))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
#The function produces a dendrogram of the eigengenes and trait(s), and a heatmap of their relationships. 
#To split the dendrogram and heatmap plots, we can use the following code
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

collectGarbage()
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")


#3.d Summary output of network analysis results
#We have found modules with high association with our trait of interest, and have identified their central players by
#the Module Membership measure. We now merge this statistical information with gene annotation and write out a
#file that summarizes the most important results and can be inspected in standard spreadsheet software such as MS
#Excel or Open Office Calc. Our expression data are only annotated by probe ID names: the command
colnames(datExprEV)
#will return all Gene IDs included in the analysis. Similarly,
colnames(datExprEV)[dynamicColorsEV=="green"]
#will return probe IDs belonging to the blue module. To facilitate interpretation of the results, 
#we will connect Uniprot IDs to gene names and universally recognized identification numbers (Entrez codes).

#will need to convert to Enterz ID
#https://yulab-smu.github.io/clusterProfiler-book/chapter14.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
#install db for mouse IDs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DOSE")
library(DOSE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ReactomePA")
library(ReactomePA)

getwd()
#create vector of IDs
View(datExprEV)
tdatExprEV<-as.data.frame(t(datExprEV))
View(tdatExprEV)
library(data.table)
setDT(tdatExprEV, keep.rownames = "Gene_ID")[]

ID <- tdatExprEV$Gene_ID
eg <- bitr(ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg)
#rename column in eg so matching title
colnames(eg)[1] <- "Gene_ID"
#merge based on Gene_ID
tdatExprEV_Entrez <- merge(eg, tdatExprEV, all.y = TRUE, by = "Gene_ID")
View(tdatExprEV_Entrez)
View(geneTraitSignificance)
View(GSPvalue)

#because 4.86% of input gene IDs are fail to map, we lost some proteins, going from 2644 to 2615 (old one was 2238 to 2124). 
#and since geneTraitSignificance and GSPvalue both have 2238, we need to merge them to a new d.f. with 2644 rows.
geneInfo0<-data.frame(Gene_ID = tdatExprEV$Gene_ID,
                      ENTREZID = tdatExprEV_Entrez$ENTREZID,
                      moduleColor = dynamicColorsEV)

View(geneInfo0)

#We now create a data frame holding the following information for all probes: probe ID, gene symbol, Locus Link ID
#(Entrez code), module color, gene significance for EV_Cell, and module membership and p-values in all modules. The
#modules will be ordered by their significance for EV_Cell, with the most significant ones to the left.
# Create the starting data frame

#previously done above:
# Define variable trait containing the Sex column of datTrait
Sex = as.data.frame(datTraitsEV$Sex);
names(Sex) = "Sex" 
# names (colors) of the modules
modNames = substring(names(MEsEV), 3)
geneModuleMembership = as.data.frame(cor(datExprEV, MEsEV, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExprEV, Sex, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Sex), sep="");
names(GSPvalue) = paste("p.GS.", names(Sex), sep="");

#up to here.

# Order modules by their significance for sample trait Sex location
modOrder = order(-abs(cor(MEsEV, Sex, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
View(geneInfo0)

#then let's get each gene's correlation and significance to the trait of choice (Sex)

geneTraitSignificance1<-geneTraitSignificance
GSPvalue1<-GSPvalue
setDT(geneTraitSignificance1, keep.rownames = "Gene_ID")[]
View(geneTraitSignificance1)
setDT(GSPvalue1, keep.rownames = "Gene_ID")[]


geneInfo0<-merge(geneInfo0,geneTraitSignificance1, by.y = "Gene_ID")
geneInfo0<-merge(geneInfo0,GSPvalue1, by.y = "Gene_ID")
#once i merge, check to see if i have different number of rows
#let's see how many unique rows i have in Protein column:
length(unique(geneInfo0$Gene_ID))
#19768 

#We now create a data frame holding the following information for all probes: probe ID, gene symbol, Locus Link ID
#(Entrez code), module color, gene significance for Sex, and module membership and p-values in all modules. The
#modules will be ordered by their significance for Sex, with the most significant ones to the left.
# Create the starting data frame



# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Sex));
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
getwd()

write.csv(geneInfo, file = "rna_WGCNA_EV_Trait_Sex_geneInfo.csv")

collectGarbage()
save.image("Transcriptomic_WCGNA.RData")

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute"))
BiocManager::install("WGCNA")
install.packages("BiocManager")
library(WGCNA)
library(Seurat)
library(flashClust)
library(cluster)
library(DOSE)
library(gplots)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#4 Interfacing network analysis with other data such as functional annotation and gene ontology
#Our previous analysis has identified several modules (labeled brown, red, and salmon) that are highly associated with weight. 
#To facilitate a biological interpretation, we would like to know the gene ontologies of the genes in the modules, 
#whether they are significantly enriched in certain functional categories etc.

#4.a Output gene lists for use with online software and services

#One option is to simply export a list of gene identifiers that can be used as input for several popular gene ontology
#and functional enrichment analysis suites such as David or AmiGO. 
#For example, we write out the LocusLinkID (entrez) codes for the brown module into a file:

View(datExprEV)
# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExprEV)
#probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = geneInfo$ENTREZID;
# $ Choose interesting modules
intModules = c("blue", "brown", "yellow", "green", "turquoise")
for (module in intModules)
{
  # Select module probes
  modGenes = (dynamicColorsEV==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
#probes = names(datExprEV)
#probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
#allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
#intModules = c("brown", "red", "salmon")
#for (module in intModules)
#{
# Select module probes
#  modGenes = (moduleColors==module)
# Get their entrez ID codes
#  modLLIDs = allLLIDs[modGenes];
# Write them into a file
#  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
#  write.table(as.data.frame(modLLIDs), file = fileName,
#              row.names = FALSE, col.names = FALSE)
#}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("EVgenes-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

GOenr = GOenrichmentAnalysis(dynamicColorsEV, allLLIDs, organism = "mouse", nBestP = 10);

tab = GOenr$bestPTerms[[4]]$enrichment

names(tab)

write.table(tab, file = "rna_WGCNA_EV_SexTrait_GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab
View(screenTab)
View(tab)
write.table(screenTab, file = "rna_WGCNA_EV_SexTrait_GOEnrichmentTable_shortened.csv", sep = ",", quote = TRUE, row.names = FALSE)

#save(MEs, moduleLabels, moduleColors, geneTree,
#     file = "FemaleLiver-02-networkConstruction-auto.RData")

save(MEs, dynamicModsEV, dynamicColorsEV, geneTreeEV,
     file = "EV-02-networkConstruction-auto.RData")


getwd()
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")


#let's try merged cut
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf



###let's use merged tree cut instead of dynamic tree cut####
#2.b.5 Merging of modules whose expression profiles are very similar
#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge
#such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we
#calculate their eigengenes and cluster them on their correlation:

#III. Using simulated data to evaluate different module detection methods
#and gene screening approaches
#6. Relating modules and module eigengenes to external data
#6.a Representing modules by eigengenes and relating eigengenes to one another
#To get a sense of how related the modules are one can summarize each module by its eigengene (first principal component).


datMEEV=moduleEigengenes(datExprEV,colorh1_EV)$eigengenes
signif(cor(datMEEV, use="p"), 2)
#The result is



#We define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation
#between the module eigengenes, and use it to cluster the eigengene:
dissimMEEV=(1-t(cor(datMEEV, method="p")))/2
hclustdatMEEV=hclust(as.dist(dissimMEEV), method="average" )
# Plot the eigengene dendrogram
sizeGrWindow(8,9)
par(mfrow=c(1,1))
plot(hclustdatMEEV, main="EV Clustering tree based of the module eigengenes (first PC)")

set.seed(082422)
#We choose a height cut of 0.2, corresponding to correlation of 0.8, to merge (see Fig. 4):
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
## Call an automatic merging function
mergeEV = mergeCloseModules(datExprEV, dynamicModsEV, cutHeight = MEDissThres, verbose = 3)
##The variable merge contains various information; we will need the following:
# Numeric module labels
moduleLabelsEV = mergeEV$colors;
# Convert labels to colors
moduleColorsEV = labels2colors(moduleLabelsEV)
# Eigengenes of the new merged modules:
mergedMEsEV = mergeEV$newMEs;

#previously:
# The merged module colors
#mergedColorsEV = mergeEV$colors;
# Eigengenes of the new merged modules:
#mergedMEsEV = mergeEV$newMEs;


#To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged
#module colors underneath (Figure 5).
sizeGrWindow(12, 9)
#pdf(file = "Plots_EV_geneDendro.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTreeEV, cbind(dynamicColorsEV, moduleColorsEV),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
table(dynamicColorsEV)
table(moduleColorsEV)
table(moduleLabelsEV)
#table(moduleColorsEV)
#black         blue        brown    darkgreen      darkred         grey      magenta midnightblue         pink       purple 
#1081         3641         4069           47         1001           17          720          295          905          622 
#salmon    turquoise       yellow 
#452         5430         1488 

#dynamicColors=staticColors
#set the diagonal of the dissimilarity to NA 
diag(dissTOMEV) = NA;
diag(dissTOMCtrlEV) = NA;
diag(dissTOM120EV) = NA;
diag(dissTOM320EV) = NA;

#change color so darker color means more connections.
library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
#dynamicColors=staticColors
#set the diagonal of the dissimilarity to NA 
#diag(dissTOM) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to a power  to bring out the module structure
#8.b Topological overlap matrix plot for visualizing the network
#We now create a so-called TOM plot, a heatmap plot depicting the topological overlap matrix 
#supplemented by hierarchical clustering dendrograms and the module colors.
#We have set the diagonal of the dissimilarity to NA and raised it to the power of 4 
#to bring out the module structure (these changes effectively amount to a change in the color scale of the plot).

#Visualize the Tom plot. Raise the dissimilarity matrix to a power  to bring out the module structure
sizeGrWindow(7,7)
#below gives you the figure in Nihal's slide
#TOMplot(dissTOMEV^4, geneTreeEV, as.character(dynamicColorsEV), main = "EV TOM Heatmap Plot, Module RNAs", cex.main=2, col=myheatcol)
#TOMplot(dissTOMCtrlEV^4, geneTreeCtrlEV, as.character(dynamicColorsCtrlEV), main = "EV_Ctrl TOM Heatmap Plot, Module RNAs", cex.main=2, col=myheatcol)
#TOMplot(dissTOM120EV^4, geneTree120EV, as.character(dynamicColors120EV), main = "EV_120 TOM Heatmap Plot, Module RNAs", col=myheatcol)
#TOMplot(dissTOM320EV^4, geneTree320EV, as.character(dynamicColors320EV), main = "EV_320 TOM Heatmap Plot, Module RNAs", col=myheatcol)

png("rna_WGCNA_EV_TOM Heatmap Plot_mergedColors_h0.2, Module RNA_sft_pwr_9.png", width = 2500, height = 2500, res = 300)
TOMplot(dissTOMEV^4, geneTreeEV, as.character(moduleColorsEV), main = "EV and Cell TOM Heatmap Plot, Module RNAs", col=myheatcol )
dev.off()

png("rna_WGCNA_EV_TOM Heatmap Plot, Module RNA_sft_pwr_9.png", width = 2500, height = 2500, res = 300)
TOMplot(dissTOMEV^4, geneTreeEV, as.character(dynamicColorsEV), main = "EV TOM Heatmap Plot, Module RNAs", col=myheatcol)
dev.off()

png("rna_WGCNA_EV_Ctrl_TOM Heatmap Plot, Module RNA_sft_pwr_18.png", width = 2500, height = 2500, res = 300)
TOMplot(dissTOMCtrlEV^4, geneTreeCtrlEV, as.character(dynamicColorsCtrlEV), main = "EV Ctrl TOM Heatmap Plot, Module RNAs", col=myheatcol)
dev.off()
png("rna_WGCNA_EV_120_TOM Heatmap Plot, Module RNA_sft_pwr_18.png", width = 2500, height = 2500, res = 300)
TOMplot(dissTOM120EV^4, geneTree120EV, as.character(dynamicColors120EV), main = "EV 120 TOM Heatmap Plot, Module RNAs", col=myheatcol)
dev.off()
png("rna_WGCNA_EV_320_TOM Heatmap Plot, Module RNA_sft_pwr_18.png", width = 2500, height = 2500, res = 300)
TOMplot(dissTOM320EV^4, geneTree320EV, as.character(dynamicColors320EV), main = "EV 320 TOM Heatmap Plot, Module RNAs", col=myheatcol)
dev.off()

#dynamicColors=staticColors
#below will give you module tables with gene names in each module
module_colorsEV= setdiff(unique(moduleColorsEV), "grey")
#module_colorsEV= setdiff(unique(dynamicColorsEV), "grey")
module_colorsCtrlEV= setdiff(unique(dynamicColorsCtrlEV), "grey")
module_colors120EV= setdiff(unique(dynamicColors120EV), "grey")
module_colors320EV= setdiff(unique(dynamicColors320EV), "grey")

#SubGeneNames=colnames(datExpr)
#SubGeneNamesEV=colnames(datExprEV)
#SubGeneNamesCtrlEV=colnames(datExprCtrlEV)
#SubGeneNames120EV=colnames(datExpr120EV)
#SubGeneNames320EV=colnames(datExpr320EV)


for (color in module_colorsEV){
  moduleEV=SubGeneNamesEV[which(moduleColorsEV==color)]
  write.table(moduleEV, paste("moduleEV_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

#Dae, run colorh1= dynamicColors
mergedcolorh1_EV= moduleColorsEV



#Dae, to get hub genes, run hubs    = chooseTopHubInEachModule(datExpr, colorh1)
mergedhubsEV    = chooseTopHubInEachModule(datExprEV, mergedcolorh1_EV)
mergedhubsEV
#       black         blue        brown    darkgreen      darkred      magenta midnightblue         pink       purple       salmon 
#"Adam18"      "Prmt1"    "Col24a1"       "Insc"        "Bik"      "Rtkn2"    "Gm34256"      "Kcnb2"       "Cps1"   "Arhgap15" 
#turquoise       yellow 
#"Cenpa"       "Sfpq" 

hubsEV    = chooseTopHubInEachModule(datExprEV, colorh1_EV)
hubsEV


###modules and traits####
getwd()
setwd("/scratch/user/daehyukchung")
#let's do a multidimensional scaling plot
#Multidimensional scaling (MDS) is a multivariate data analysis approach 
#that is used to visualize the similarity/dissimilarity between samples by plotting points in two dimensional plots.
cmd1EV=cmdscale(as.dist(dissTOMEV),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd1EV, col=as.character(mergedcolorh1_EV), main="EV MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")


######




####modules to external traits relationship###
#here, we want to see which modules relate closely to which traits
#for my samples, traits will be pregnancy/set, sex, and ethanol.
#it can also be ev or cell, if i have all the samples together.
#Module-Trait relationships. Color scale (red-green) represents the strength of the correlation between the module and the trait. 
#Each box gives a correlation value (R^2) followed by p-value (in parenthesis). 
#the color scale is the strength of the correlation between the module and the trait, 
#where for example, MEblack and Female being very red means highly correlated.

#let's load the excel csv file i have created for traits. 
#i am adding "row.names=1" function, which means make the first column (sample names) into row names 
#datTraits = read.csv("Traits_wgcna1.csv",row.names = 1)
#write.csv(datTraits,file="Traitswgcna.csv")
#form a data frame analogous to expression data that will hold the clinical traits.
#let's match the rownames of datExpr and datTraits, after you check that the samples are lined up correctly by rows
rownames(datTraits)
#let's get row 1:18 that only have EV samples, and column 2:4 to exclude EV_Cell trait
datTraitsEV<-datTraits[1:18,2:4]
rownames(datTraitsEV)
rownames(datExprEV)


rownames(datExprEV)<-rownames(datTraitsEV)

table(rownames(datTraitsEV)==rownames(datExprEV)) 
#should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraitsEV)

##
getwd()

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


#3 Relating modules to external clinical traits
#3.a Quantifying module-trait associations
#In this analysis we would like to identify modules that are significantly associated with the measured clinical traits.
#Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external
#traits and look for the most significant associations:

#instead of moduleColors, i have dynamicColorsEV from:
#dynamicColorsEV = labels2colors(dynamicModsEV)
View(dynamicColorsEV)
View(colorh1)

# Define numbers of genes and samples
nGenesEV = ncol(datExprEV);
nSamplesEV = nrow(datExprEV);
# Recalculate MEs with color labels
#moduleEigengenes() simply calculate the 1st Principal Component (PC) i.e., module eigengene (ME), of each module.
mergedMEs0EV = moduleEigengenes(datExprEV, moduleColorsEV)$eigengenes
mergedMEs0EV = orderMEs(mergedMEs0EV)
mergedmoduleTraitCorEV = cor(mergedMEs0EV, datTraitsEV, use = "p");
mergedmoduleTraitPvalueEV = corPvalueStudent(mergedmoduleTraitCorEV, nSamplesEV);


#let's use a graphical representation to read the table, a color-coded table,
#where we color code each association of module eigengenes to the sample trait by the correlation value:
sizeGrWindow(10,6)
sizeGrWindow(16,10)
# Will display correlations and their p-values
textMatrix = paste(signif(mergedmoduleTraitCorEV, 2), "\n(",
                   signif(mergedmoduleTraitPvalueEV, 1), ")", sep = "");
dim(textMatrix) = dim(mergedmoduleTraitCorEV)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = mergedmoduleTraitCorEV,
               xLabels = names(datTraitsEV),
               yLabels = names(mergedMEs0EV),
               ySymbols = names(mergedMEs0EV),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("EV Module-Trait Relationships"))

#let's match the rownames of datExpr and datTraits, after you check that the samples are lined up correctly by rows
rownames(datExprEV)<-rownames(datTraitsEV)

table(rownames(datTraitsEV)==rownames(datExprEV)) 
#should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraitsEV)
y=datTraitsEV$Sex
y

mergeddatMEEV=moduleEigengenes(datExprEV,mergedcolorh1_EV)$eigengenes
signif(cor(mergeddatMEEV, use="p"), 2)
#The result is

#6.e.2 Measure of module significance as average gene significance
#One can also define a measure of module significance as the average gene significance of all genes in the module. We
#use the absolute value for defining a correlation based gene significance measure.
GS1=as.numeric(cor(y,datExprEV, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
mergedModuleSignificanceEV=tapply(GeneSignificance, mergedcolorh1_EV, mean, na.rm=T)


mergedModuleSignificanceEV
#       black         blue        brown    darkgreen      darkred         grey      magenta midnightblue         pink       purple 
#0.2933319    0.2459879    0.3457895    0.3246322    0.1810582    0.5268715    0.2821780    0.2806292    0.3354527    0.2636081 
#salmon    turquoise       yellow 
#0.4597463    0.5130779    0.1815087 

#To plot module significance, one can use
sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,mergedcolorh1_EV)

#We define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation
#between the module eigengenes, and use it to cluster the eigengene:
mergeddissimMEEV=(1-t(cor(mergeddatMEEV, method="p")))/2
mergedhclustdatMEEV=hclust(as.dist(mergeddissimMEEV), method="average" )
# Plot the eigengene dendrogram
sizeGrWindow(8,9)
par(mfrow=c(1,1))
plot(mergedhclustdatMEEV, main="EV Clustering tree based of the module eigengenes (first PC)")



#The advantage of this approach is that it can be used for any gene significance measure. 
#A gene significance measure could be defined without reference to a sample trait. 
#For example, it could indicate pathway membership (1 or 0) or gene essentiality (1 or 0), etc.
#The next logical step in an analysis of empirical data would be to carry out a functional enrichment analysis of
#each module, for example using the software EASE (David) at http://david.abcc.ncifcrf.gov/summary.jsp.
#We refer the reader to Tutorial I for an example of a functional enrichment analysis. Before we end, we again save
#calculated results for use in subsequent sections.
collectGarbage()
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")
getwd()

#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf

#3.b Gene relationship to trait and important modules: Gene Significance and Module Membership
#We quantify associations of individual genes with our trait of interest (Sex) by defining Gene Significance GS as
#(the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative
#measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 
#This allows us to quantify the similarity of all genes on the array to every module.
# Define variable trait containing the Sex column of datTrait
Sex = as.data.frame(datTraitsEV$Sex);
names(Sex) = "Sex"
# names (colors) of the modules
modNames = substring(names(mergedMEs0EV), 3)
geneModuleMembership = as.data.frame(cor(datExprEV, mergedMEs0EV, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExprEV, Sex, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Sex), sep="");
names(GSPvalue) = paste("p.GS.", names(Sex), sep="");


#3.c Intramodular analysis: identifying genes with high GS and MM
#Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module
#membership in interesting modules. As an example, we look at the lightcyan, magenta, purple, turquoise, blue module that has the highest association
#with Sex We plot a scatterplot of Gene Significance vs. Module Membership in the blue module:
#module lightcyon, magenta, purple highest gene significance for female, turquoise, blue for male.
# names (colors) of the modules
library(ggplot2)

module = "lightcyan"
column = match(module, modNames);
moduleGenes = moduleColorsEV==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "cyan")

module = "magenta"
column = match(module, modNames);
moduleGenes = moduleColorsEV==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "purple"
column = match(module, modNames);
moduleGenes = dynamicColorsEV==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module)


module = "turquoise"
column = match(module, modNames);
moduleGenes = dynamicColorsEV==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "blue"
column = match(module, modNames);
moduleGenes = dynamicColorsEV==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



#I. Network analysis of liver expression data in female mice
#5. Network visualization using WGCNA functions
#5.b Visualizing the network of eigengenes
#It is often interesting to study the relationships among the found modules. 
#One can use the eigengenes as representative profiles and quantify module similarity by eigengene correlation. 
#The package contains a convenient function plotEigengeneNetworks that generates a summary plot of the eigengene network. 
#It is usually informative to add a clinical trait (or multiple traits) to the eigengenes 
#to see how the traits fit into the eigengene network:
#Previously: 
#MEs0 = moduleEigengenes(datExprCell, dynamicColorsCell)$eigengenes
#MEs = orderMEs(MEs0)
## Isolate weight from the clinical traits
#Sex = as.data.frame(datTraits$Sex);
#names(Sex) = "Sex"
# Add the Sex to existing module eigengenes
MET = orderMEs(cbind(mergedMEs0EV, Sex))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
#The function produces a dendrogram of the eigengenes and trait(s), and a heatmap of their relationships. 
#To split the dendrogram and heatmap plots, we can use the following code
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

collectGarbage()
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")


#3.d Summary output of network analysis results
#We have found modules with high association with our trait of interest, and have identified their central players by
#the Module Membership measure. We now merge this statistical information with gene annotation and write out a
#file that summarizes the most important results and can be inspected in standard spreadsheet software such as MS
#Excel or Open Office Calc. Our expression data are only annotated by probe ID names: the command
colnames(datExprEV)
#will return all Gene IDs included in the analysis. Similarly,
colnames(datExprEV)[dynamicColorsEV=="green"]
#will return probe IDs belonging to the blue module. To facilitate interpretation of the results, 
#we will connect Uniprot IDs to gene names and universally recognized identification numbers (Entrez codes).

#will need to convert to Enterz ID
#https://yulab-smu.github.io/clusterProfiler-book/chapter14.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
#install db for mouse IDs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DOSE")
library(DOSE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ReactomePA")
library(ReactomePA)

getwd()
#create vector of IDs
View(datExprEV)
tdatExprEV<-as.data.frame(t(datExprEV))
View(tdatExprEV)
library(data.table)
setDT(tdatExprEV, keep.rownames = "Gene_ID")[]

ID <- tdatExprEV$Gene_ID
eg <- bitr(ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg)
#rename column in eg so matching title
colnames(eg)[1] <- "Gene_ID"
#merge based on Gene_ID
tdatExprEV_Entrez <- merge(eg, tdatExprEV, all.y = TRUE, by = "Gene_ID")
View(tdatExprEV_Entrez)
View(geneTraitSignificance)
View(GSPvalue)

#because 4.86% of input gene IDs are fail to map, we lost some proteins, going from 2644 to 2615 (old one was 2238 to 2124). 
#and since geneTraitSignificance and GSPvalue both have 2238, we need to merge them to a new d.f. with 2644 rows.
geneInfo0<-data.frame(Gene_ID = tdatExprEV$Gene_ID,
                      ENTREZID = tdatExprEV_Entrez$ENTREZID,
                      moduleColor = moduleColorsEV)

View(geneInfo0)

# Order modules by their significance for sample trait Sex location
modOrder = order(-abs(cor(mergedMEs0EV, Sex, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
View(geneInfo0)



geneTraitSignificance1<-geneTraitSignificance
GSPvalue1<-GSPvalue
setDT(geneTraitSignificance1, keep.rownames = "Gene_ID")[]
View(geneTraitSignificance1)
setDT(GSPvalue1, keep.rownames = "Gene_ID")[]



geneInfo0<-merge(geneInfo0,geneTraitSignificance1, by.y = "Gene_ID")
geneInfo0<-merge(geneInfo0,GSPvalue1, by.y = "Gene_ID")
#once i merge, check to see if i have different number of rows
#let's see how many unique rows i have in Protein column:
length(unique(geneInfo0$Gene_ID))
#19768 
#if i have duplicates,
#remove them. 
geneInfo0<-geneInfo0[!duplicated(geneInfo0$Gene_ID), ]

#let's move the columns around
geneInfo0<-geneInfo0[,c(1,2,30,31,3:29)]
View(geneInfo0)


# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Sex));
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
getwd()

write.csv(geneInfo, file = "rna_WGCNA_EV_geneInfo_sft_pwr_9_h0.2.csv")

#let's pick columns with Gene_ID, ENTREZID, and moduleColor
rna_WGCNA_EV_Gene_ID_moduleColor_sft_pwr_9_h0.2<-
  geneInfo[,c("Gene_ID","ENTREZID","moduleColor")]
write.csv(rna_WGCNA_EV_Gene_ID_moduleColor_sft_pwr_9_h0.2, file = "rna_WGCNA_EV_Gene_ID_moduleColor_sft_pwr_9_h0.2.csv")

collectGarbage()

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute"))
BiocManager::install("WGCNA")
install.packages("BiocManager")
library(WGCNA)
library(Seurat)
library(flashClust)
library(cluster)
library(DOSE)
library(gplots)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#4 Interfacing network analysis with other data such as functional annotation and gene ontology
#Our previous analysis has identified several modules (labeled brown, red, and salmon) that are highly associated with weight. 
#To facilitate a biological interpretation, we would like to know the gene ontologies of the genes in the modules, 
#whether they are significantly enriched in certain functional categories etc.

#4.a Output gene lists for use with online software and services

#One option is to simply export a list of gene identifiers that can be used as input for several popular gene ontology
#and functional enrichment analysis suites such as David or AmiGO. 
#For example, we write out the LocusLinkID (entrez) codes for the brown module into a file:

View(datExprEV)
# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExprEV)
#probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = geneInfo$ENTREZID;
# $ Choose interesting modules
intModules = c("blue", "brown", "yellow", "green", "turquoise")
for (module in intModules)
{
  # Select module probes
  modGenes = (dynamicColorsEV==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
#probes = names(datExprEV)
#probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
#allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
#intModules = c("brown", "red", "salmon")
#for (module in intModules)
#{
# Select module probes
#  modGenes = (moduleColors==module)
# Get their entrez ID codes
#  modLLIDs = allLLIDs[modGenes];
# Write them into a file
#  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
#  write.table(as.data.frame(modLLIDs), file = fileName,
#              row.names = FALSE, col.names = FALSE)
#}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("EVgenes-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

GOenr = GOenrichmentAnalysis(moduleColorsEV, allLLIDs, organism = "mouse", nBestP = 10);

tab = GOenr$bestPTerms[[4]]$enrichment

names(tab)

write.table(tab, file = "rna_WGCNA_EV_SexTrait_GOEnrichmentTable_sft_pwr_9_h0.2.csv", sep = ",", quote = TRUE, row.names = FALSE)

keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab
View(screenTab)
View(tab)
write.table(screenTab, file = "rna_WGCNA_EV_SexTrait_GOEnrichmentTable_shortened_sft_pwr_9_h0.2.csv", sep = ",", quote = TRUE, row.names = FALSE)


#previously:
#save(MEs, moduleLabels, moduleColors, geneTree,
#     file = "FemaleLiver-02-networkConstruction-auto.RData")
##be performed automatically:
#merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)  
#The variable merge contains various information; we will need the following:
# Numeric module labels
#moduleLabels = merge$colors;
# Convert labels to colors
#moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
#consMEs = merge$newMEs;
#
## Call an automatic merging function
#mergeEV = mergeCloseModules(datExprEV, dynamicColorsEV, cutHeight = MEDissThres, verbose = 3)
##The variable merge contains various information; we will need the following:
# Numeric module labels
#moduleLabelsEV = mergeEV$colors;
# Convert labels to colors
#moduleColorsEV = labels2colors(moduleLabelsEV)
# Eigengenes of the new merged modules:
#mergedMEsEV = mergeEV$newMEs;

table(dynamicColorsEV)
table(moduleColors)
table(dynamicLabels)
table(moduleLabels)
table(mergeEV$colors)
table(labels2colors(moduleLabelsEV))
#
save(mergedMEsEV, moduleLabelsEV, moduleColorsEV, geneTreeEV,
     file = "EV-02-networkConstruction-man.RData")

#save(MEs, dynamicModsEV, dynamicColorsEV, geneTreeEV,
#     file = "EV-02-networkConstruction-auto.RData")


getwd()
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")

######





########nice. now, lets do it for cells.####### 
library(WGCNA)
library(Seurat)
library(flashClust)
library(cluster)
library(DOSE)
library(gplots)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()


datExprCell<-datExpr[19:36,]


##which column has gene sum that has 0 value 

which(colSums(datExprCell, na.rm = TRUE)==0)


#5.a Defining a weighted gene co-expression network
# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1Cell=abs(cor(datExprCell,use="p"))^6

# When you have relatively few genes (<5000) use the following code. 
#k is for connectivity
#kCell=as.vector(apply(ADJ1Cell,2,sum, na.rm=T))


# When you have a lot of genes use the following code
#datE<-datExpr
#k=softConnectivity(datE,power=6)
kCell=softConnectivity(datExprCell,power=6)

# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(kCell)
scaleFreePlot(kCell, main="Check scale free topology\n")
hist(kCtrlCell)

#5.c.2 Use of topologial overlap to define dissimilarity
#Adjacency can be used to define a separate measure of similarity, the Topological Overlap Matrix(TOM) [2, 1]:
dissTOMCell=TOMdist(ADJ1Cell)

collectGarbage()

getwd()
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")

#let's get gene names from our data

SubGeneNamesCell=colnames(datExprCell)


getwd()
setwd("/scratch/user/daehyukchung")
#choose softpower
powers = c(c(1:10), seq(from = 12, to=20, by=2));
softPowerCell=pickSoftThreshold(datExprCell,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "signed")

#previous issue that has now been reserved by excluding certain proteins before making datExpr:
##issues with ctrl cell samples, where i Have a 0 value for the entire column.
#is.na(datExprCtrlCell[,1])
#length(is.na(datExprCtrlCell[,1]))
##which column has gene sum that has 0 value 

which(colSums(datExprCell, na.rm = TRUE)==0)


#
#Plot the results
sftCell<-softPowerCell


#to get 2 separate plots in one window:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sftCell$fitIndices[,1], -sign(sftCell$fitIndices[,3])*sftCell$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sftCell$fitIndices[,1], -sign(sftCell$fitIndices[,3])*sftCell$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sftCell$fitIndices[,1], sftCell$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sftCell$fitIndices[,1], sftCell$fitIndices[,5], labels=powers, cex=cex1,col="red")

#so, looking at the mean connectivity, they are not all similar, 
#you usually want to choose a soft threshold (power) that is over 0.8 for R^2,
#which would be 14 for 36 EV and cell samples, so lets try that.
#for EV, it is 6 and for Cell, it is 9. The reason 36 samples of EV and Cell is 14 is due to the fact that 
#EV and Cell sample groups are so different from each other. 
#so use softPower 14 for the 36 samples, and use softPower 9 for EV and Cell so that we can compare EV with Cell. 
#the one above 0.8 line is 6, 14, none (too high), 20.
#the appropriate soft-thresholding power can be chosen based on the number of samples as in the table below. 
#This table has been updated in December 2017 to make the resulting networks conservative.
#Number of samples	Unsigned and signed hybrid networks	Signed networks
#Less than 20	                  9	                            18
#20-30	                        8	                            16
#30-40	                        7	                            14
#more than 40	                  6	                            12
#for more info: https://bioinformatics.stackexchange.com/questions/11335/wgcna-problem-with-selecting-soft-threshold

softPowerCell = 9



#let's do a adjacency matrix (kind of like a correlation matrix)
#adj= adjacency(datExprCell,type = "signed", power = softPower);

#turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
#trying to calculate how connected genes are to each other.
TOMCell=TOMsimilarityFromExpr(datExprCell,networkType = "signed", TOMType = "signed", power = softPowerCell);


#since correlation matrix will not have names, we are giving subgenenames to column and row names. 
#because we removed a column named Q8BYM5 for ctrl group, we do:


colnames(TOMCell) =rownames(TOMCell) =SubGeneNamesCell
dissTOMCell=1-TOMCell


#
#hierarchical clustering of the genes based on the TOM dissimilarity measure
library(flashClust)
geneTreeCell = flashClust(as.dist(dissTOMCell),method="average")

#plot the resulting clustering tree (dendrogram)
par(mfrow = c(1,1))
plot(geneTreeCell, xlab="", sub="",cex=0.3,main=paste0("genesCell"))
#from that, we see how the genes are clustered together. 
#we will cut the tree at a horizontal line, let's say 0.8 height, then whatever is cut below will form a module.
#dynamiccuttree will do this automatically for you.

# Set the minimum module size
minModuleSize = 30;

# Module identification using dynamic tree cut
#important to identify modules
dynamicModsCell = cutreeDynamic(dendro = geneTreeCell,  method="tree", minClusterSize = minModuleSize);

# Plot the dendrogram with module colors

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
#tells how many genes i have in each module and how many modules there are. 
table(dynamicModsCtrlCell)

#table(colorStaticADJ)
##Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColorsCell = labels2colors(dynamicModsCell)
table(dynamicColorsCell)

#staticColors = labels2colors(colorStaticADJ)
#table(staticColors)
#
plotDendroAndColors(geneTreeCell, dynamicColorsCell, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Cell Gene dendrogram and module colors")
#so turquoise color is leftover genes that dont really belong to any module
#the other modules are important, which we will extract and study more later. 
#


###let's use merged tree cut instead of dynamic tree cut####
#2.b.5 Merging of modules whose expression profiles are very similar
#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge
#such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we
#calculate their eigengenes and cluster them on their correlation:

#III. Using simulated data to evaluate different module detection methods
#and gene screening approaches
#6. Relating modules and module eigengenes to external data
#6.a Representing modules by eigengenes and relating eigengenes to one another
#To get a sense of how related the modules are one can summarize each module by its eigengene (first principal component).
colorh1_Cell= dynamicColorsCell

datMECell=moduleEigengenes(datExprCell,colorh1_Cell)$eigengenes
signif(cor(datMECell, use="p"), 2)
#The result is


#We define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation
#between the module eigengenes, and use it to cluster the eigengene:
dissimMECell=(1-t(cor(datMECell, method="p")))/2
hclustdatMECell=hclust(as.dist(dissimMECell), method="average" )
# Plot the eigengene dendrogram
sizeGrWindow(8,9)
par(mfrow=c(1,1))
plot(hclustdatMECell, main="Cell Clustering tree based of the module eigengenes (first PC)")

set.seed(082522)
#We choose a height cut of 0.2, corresponding to correlation of 0.8, to merge (see Fig. 4):
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
## Call an automatic merging function
mergeCell = mergeCloseModules(datExprCell, dynamicModsCell, cutHeight = MEDissThres, verbose = 3)
##The variable merge contains various information; we will need the following:
# Numeric module labels
moduleLabelsCell = mergeCell$colors;
# Convert labels to colors
moduleColorsCell = labels2colors(moduleLabelsCell)
# Eigengenes of the new merged modules:
mergedMEsCell = mergeCell$newMEs;

#previously:
# The merged module colors
#mergedColorsEV = mergeEV$colors;
# Eigengenes of the new merged modules:
#mergedMEsEV = mergeEV$newMEs;


#To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged
#module colors underneath (Figure 5).
sizeGrWindow(12, 9)
#pdf(file = "Plots_Cell_geneDendro.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTreeCell, cbind(dynamicColorsCell, moduleColorsCell),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
table(dynamicColorsCell)
table(moduleColorsCell)
table(moduleLabelsCell)
#table(moduleColorsCell)
#black        blue       brown   darkgreen       green        grey   lightcyan 
#2776        2889        2108          83        1675          51         249 
#lightyellow     magenta         red      yellow 
#717        5983        1383        1854 

#dynamicColors=staticColors
#set the diagonal of the dissimilarity to NA 
diag(dissTOMCell) = NA;


#change color so darker color means more connections.
library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
#dynamicColors=staticColors
#set the diagonal of the dissimilarity to NA 
#diag(dissTOM) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to a power  to bring out the module structure
#8.b Topological overlap matrix plot for visualizing the network
#We now create a so-called TOM plot, a heatmap plot depicting the topological overlap matrix 
#supplemented by hierarchical clustering dendrograms and the module colors.
#We have set the diagonal of the dissimilarity to NA and raised it to the power of 4 
#to bring out the module structure (these changes effectively amount to a change in the color scale of the plot).

#Visualize the Tom plot. Raise the dissimilarity matrix to a power  to bring out the module structure

png("rna_WGCNA_Cell_TOM Heatmap Plot_mergedColors_h0.2, Module RNA_sft_pwr_9.png", width = 2500, height = 2500, res = 300)
TOMplot(dissTOMCell^4, geneTreeCell, as.character(moduleColorsCell), main = "Cell TOM Heatmap Plot, Module RNAs", col=myheatcol )
dev.off()

png("rna_WGCNA_Cell_TOM Heatmap Plot, Module RNA_sft_pwr_9.png", width = 2500, height = 2500, res = 300)
TOMplot(dissTOMCell^4, geneTreeCell, as.character(dynamicColorsCell), main = "Cell TOM Heatmap Plot, Module RNAs", col=myheatcol)
dev.off()



#dynamicColors=staticColors
#below will give you module tables with gene names in each module
module_colorsCell= setdiff(unique(moduleColorsCell), "grey")


for (color in module_colorsCell){
  moduleCell=SubGeneNamesCell[which(moduleColorsCell==color)]
  write.table(moduleCell, paste("moduleCell_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}


mergedcolorh1_Cell= moduleColorsCell

colorh1_Cell= dynamicColorsCell


#Dae, to get hub genes, run hubs    = chooseTopHubInEachModule(datExpr, colorh1)
mergedhubsCell    = chooseTopHubInEachModule(datExprCell, mergedcolorh1_Cell)
mergedhubsCell
#      black        blue       brown   darkgreen       green   lightcyan lightyellow 
#"Smarca2"    "Fbxo34"      "Crcp"      "Phyh"  "n-TTagt6"      "Sgcb"     "Actn4" 
#magenta         red      yellow 
#"Ssrp1"      "Lhx4"    "Parp11" 

hubsCell    = chooseTopHubInEachModule(datExprCell, colorh1_Cell)
hubsCell


#mergedhubsCell
#black         blue        brown    darkgreen      darkred      magenta midnightblue         pink       purple       salmon 
#"Adam18"      "Prmt1"    "Col24a1"       "Insc"        "Bik"      "Rtkn2"    "Gm34256"      "Kcnb2"       "Cps1"   "Arhgap15" 
#turquoise       yellow 
#"Cenpa"       "Sfpq" 

#> hubsCell
#black            blue           brown            cyan       darkgreen        darkgrey      darkorange 
#"Ttc29"         "Prmt1"        "Hs6st3"    "D8Ertd738e"          "Dpyd"        "Fam20a"          "Insc" 
#darkred   darkturquoise           green     greenyellow          grey60       lightcyan      lightgreen 
#"Bik"         "Prkcq"        "Anapc4"       "Nectin4"         "Atp5e" "1810007D17Rik"        "Lrrtm4" 
#lightyellow         magenta    midnightblue          orange            pink          purple             red 
#"Lhfpl3"         "Rtkn2"      "Serpini1"        "Gm1043"        "Ccser1"          "Cps1"       "Ccdc180" 
#royalblue          salmon         skyblue             tan       turquoise           white          yellow 
#"Mroh2b"       "Gm34256"       "Gm20471"        "Ifnar2"          "Scd2"         "Foxf2"          "Sfpq" 


###modules and traits####
getwd()
setwd("/scratch/user/daehyukchung")
#let's do a multidimensional scaling plot
#Multidimensional scaling (MDS) is a multivariate data analysis approach 
#that is used to visualize the similarity/dissimilarity between samples by plotting points in two dimensional plots.
cmd1Cell=cmdscale(as.dist(dissTOMCell),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd1Cell, col=as.character(mergedcolorh1_Cell), main="Cell MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")


######




####modules to external traits relationship###
#here, we want to see which modules relate closely to which traits
#for my samples, traits will be pregnancy/set, sex, and ethanol.
#it can also be ev or cell, if i have all the samples together.
#Module-Trait relationships. Color scale (red-green) represents the strength of the correlation between the module and the trait. 
#Each box gives a correlation value (R^2) followed by p-value (in parenthesis). 
#the color scale is the strength of the correlation between the module and the trait, 
#where for example, MEblack and Female being very red means highly correlated.

#let's load the excel csv file i have created for traits. 
#i am adding "row.names=1" function, which means make the first column (sample names) into row names 
#datTraits = read.csv("Traits_wgcna1.csv",row.names = 1)
#write.csv(datTraits,file="Traitswgcna.csv")
#form a data frame analogous to expression data that will hold the clinical traits.
#let's match the rownames of datExpr and datTraits, after you check that the samples are lined up correctly by rows
rownames(datTraits)
#let's get row 1:18 that only have Cell samples, and column 2:4 to exclude EV_Cell trait
datTraitsCell<-datTraits[19:36,2:4]
rownames(datTraitsCell)
rownames(datExprCell)


rownames(datExprCell)<-rownames(datTraitsCell)

table(rownames(datTraitsCell)==rownames(datExprCell)) 
#should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraitsCell)

##
getwd()

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


#3 Relating modules to external clinical traits
#3.a Quantifying module-trait associations
#In this analysis we would like to identify modules that are significantly associated with the measured clinical traits.
#Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external
#traits and look for the most significant associations:

#instead of moduleColors, i have dynamicColorsCell from:
#dynamicColorsCell = labels2colors(dynamicModsCell)
View(dynamicColorsCell)
View(colorh1)

# Define numbers of genes and samples
nGenesCell = ncol(datExprCell);
nSamplesCell = nrow(datExprCell);
# Recalculate MEs with color labels
#moduleEigengenes() simply calculate the 1st Principal Component (PC) i.e., module eigengene (ME), of each module.
mergedMEs0Cell = moduleEigengenes(datExprCell, moduleColorsCell)$eigengenes
mergedMEs0Cell = orderMEs(mergedMEs0Cell)
mergedmoduleTraitCorCell = cor(mergedMEs0Cell, datTraitsCell, use = "p");
mergedmoduleTraitPvalueCell = corPvalueStudent(mergedmoduleTraitCorCell, nSamplesCell);


#let's use a graphical representation to read the table, a color-coded table,
#where we color code each association of module eigengenes to the sample trait by the correlation value:
sizeGrWindow(10,6)
sizeGrWindow(16,10)
# Will display correlations and their p-values
textMatrix = paste(signif(mergedmoduleTraitCorCell, 2), "\n(",
                   signif(mergedmoduleTraitPvalueCell, 1), ")", sep = "");
dim(textMatrix) = dim(mergedmoduleTraitCorCell)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = mergedmoduleTraitCorCell,
               xLabels = names(datTraitsCell),
               yLabels = names(mergedMEs0Cell),
               ySymbols = names(mergedMEs0Cell),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Cell Module-Trait Relationships"))

#let's match the rownames of datExpr and datTraits, after you check that the samples are lined up correctly by rows
rownames(datExprCell)<-rownames(datTraitsCell)

table(rownames(datTraitsCell)==rownames(datExprCell)) 
#should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraitsCell)
y=datTraitsCell$Sex
y

mergeddatMECell=moduleEigengenes(datExprCell,mergedcolorh1_Cell)$eigengenes
signif(cor(mergeddatMECell, use="p"), 2)
#The result is

#6.e.2 Measure of module significance as average gene significance
#One can also define a measure of module significance as the average gene significance of all genes in the module. We
#use the absolute value for defining a correlation based gene significance measure.
GS1=as.numeric(cor(y,datExprCell, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
mergedModuleSignificanceCell=tapply(GeneSignificance, mergedcolorh1_Cell, mean, na.rm=T)


mergedModuleSignificanceCell
# black        blue       brown   darkgreen       green        grey   lightcyan 
#0.2186194   0.4972008   0.2553563   0.1645306   0.3624965   0.1991443   0.4546617 
#lightyellow     magenta         red      yellow 
#0.5359865   0.1938239   0.5433448   0.2907327

#To plot module significance, one can use
sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,mergedcolorh1_Cell)

#We define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation
#between the module eigengenes, and use it to cluster the eigengene:
mergeddissimMECell=(1-t(cor(mergeddatMECell, method="p")))/2
mergedhclustdatMECell=hclust(as.dist(mergeddissimMECell), method="average" )
# Plot the eigengene dendrogram
sizeGrWindow(8,9)
par(mfrow=c(1,1))
plot(mergedhclustdatMECell, main="Cell Clustering tree based of the module eigengenes (first PC)")



#The advantage of this approach is that it can be used for any gene significance measure. 
#A gene significance measure could be defined without reference to a sample trait. 
#For example, it could indicate pathway membership (1 or 0) or gene essentiality (1 or 0), etc.
#The next logical step in an analysis of empirical data would be to carry out a functional enrichment analysis of
#each module, for example using the software EASE (David) at http://david.abcc.ncifcrf.gov/summary.jsp.
#We refer the reader to Tutorial I for an example of a functional enrichment analysis. Before we end, we again save
#calculated results for use in subsequent sections.
collectGarbage()
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")
getwd()

#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf

#3.b Gene relationship to trait and important modules: Gene Significance and Module Membership
#We quantify associations of individual genes with our trait of interest (Sex) by defining Gene Significance GS as
#(the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative
#measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 
#This allows us to quantify the similarity of all genes on the array to every module.
# Define variable trait containing the Sex column of datTrait
Sex = as.data.frame(datTraitsCell$Sex);
names(Sex) = "Sex"
# names (colors) of the modules
modNames = substring(names(mergedMEs0Cell), 3)
geneModuleMembership = as.data.frame(cor(datExprCell, mergedMEs0Cell, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExprCell, Sex, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Sex), sep="");
names(GSPvalue) = paste("p.GS.", names(Sex), sep="");


#3.c Intramodular analysis: identifying genes with high GS and MM
#Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module
#membership in interesting modules. As an example, we look at the lightcyan, magenta, purple, turquoise, blue module that has the highest association
#with Sex We plot a scatterplot of Gene Significance vs. Module Membership in the blue module:
#module lightcyon, magenta, purple highest gene significance for female, turquoise, blue for male.
# names (colors) of the modules
library(ggplot2)

module = "lightcyan"
column = match(module, modNames);
moduleGenes = moduleColorsCell==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "cyan")

module = "magenta"
column = match(module, modNames);
moduleGenes = moduleColorsCell==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "purple"
column = match(module, modNames);
moduleGenes = dynamicColorsCell==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col =module)


module = "turquoise"
column = match(module, modNames);
moduleGenes = dynamicColorsCell==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "blue"
column = match(module, modNames);
moduleGenes = dynamicColorsCell==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Sex location",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



#I. Network analysis of liver expression data in female mice
#5. Network visualization using WGCNA functions
#5.b Visualizing the network of eigengenes
#It is often interesting to study the relationships among the found modules. 
#One can use the eigengenes as representative profiles and quantify module similarity by eigengene correlation. 
#The package contains a convenient function plotEigengeneNetworks that generates a summary plot of the eigengene network. 
#It is usually informative to add a clinical trait (or multiple traits) to the eigengenes 
#to see how the traits fit into the eigengene network:
#Previously: 
#MEs0 = moduleEigengenes(datExprCell, dynamicColorsCell)$eigengenes
#MEs = orderMEs(MEs0)
## Isolate weight from the clinical traits
#Sex = as.data.frame(datTraits$Sex);
#names(Sex) = "Sex"
# Add the Sex to existing module eigengenes
MET = orderMEs(cbind(mergedMEs0Cell, Sex))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
#The function produces a dendrogram of the eigengenes and trait(s), and a heatmap of their relationships. 
#To split the dendrogram and heatmap plots, we can use the following code
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

collectGarbage()
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")


#3.d Summary output of network analysis results
#We have found modules with high association with our trait of interest, and have identified their central players by
#the Module Membership measure. We now merge this statistical information with gene annotation and write out a
#file that summarizes the most important results and can be inspected in standard spreadsheet software such as MS
#Excel or Open Office Calc. Our expression data are only annotated by probe ID names: the command
colnames(datExprCell)
#will return all Gene IDs included in the analysis. Similarly,
colnames(datExprCell)[dynamicColorsCell=="green"]
#will return probe IDs belonging to the blue module. To facilitate interpretation of the results, 
#we will connect Uniprot IDs to gene names and universally recognized identification numbers (Entrez codes).

#will need to convert to Enterz ID
#https://yulab-smu.github.io/clusterProfiler-book/chapter14.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
#install db for mouse IDs
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DOSE")
library(DOSE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ReactomePA")
library(ReactomePA)

getwd()
#create vector of IDs
View(datExprCell)
tdatExprCell<-as.data.frame(t(datExprCell))
View(tdatExprCell)
library(data.table)
setDT(tdatExprCell, keep.rownames = "Gene_ID")[]

ID <- tdatExprCell$Gene_ID
eg <- bitr(ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(eg)
#rename column in eg so matching title
colnames(eg)[1] <- "Gene_ID"
#merge based on Gene_ID
tdatExprCell_Entrez <- merge(eg, tdatExprCell, all.y = TRUE, by = "Gene_ID")
View(tdatExprCell_Entrez)
View(geneTraitSignificance)
View(GSPvalue)

#because 4.86% of input gene IDs are fail to map, we lost some proteins, going from 2644 to 2615 (old one was 2238 to 2124). 
#and since geneTraitSignificance and GSPvalue both have 2238, we need to merge them to a new d.f. with 2644 rows.
geneInfo0<-data.frame(Gene_ID = tdatExprCell$Gene_ID,
                      ENTREZID = tdatExprCell_Entrez$ENTREZID,
                      moduleColor = moduleColorsCell)

View(geneInfo0)



#We now create a data frame holding the following information for all probes: probe ID, gene symbol, Locus Link ID
#(Entrez code), module color, gene significance for Sex, and module membership and p-values in all modules. The
#modules will be ordered by their significance for Sex, with the most significant ones to the left.
# Create the starting data frame

# Order modules by their significance for sample trait Sex location
modOrder = order(-abs(cor(mergedMEs0Cell, Sex, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
View(geneInfo0)

geneTraitSignificance1<-geneTraitSignificance
GSPvalue1<-GSPvalue
setDT(geneTraitSignificance1, keep.rownames = "Gene_ID")[]
View(geneTraitSignificance1)
setDT(GSPvalue1, keep.rownames = "Gene_ID")[]



geneInfo0<-merge(geneInfo0,geneTraitSignificance1, by.y = "Gene_ID")
geneInfo0<-merge(geneInfo0,GSPvalue1, by.y = "Gene_ID")
#once i merge, check to see if i have different number of rows
#let's see how many unique rows i have in Protein column:
length(unique(geneInfo0$Gene_ID))
#19768 
#if i have duplicates,
#remove them. 
geneInfo0<-geneInfo0[!duplicated(geneInfo0$Gene_ID), ]
#let's move the columns around
geneInfo0<-geneInfo0[,c(1,2,26,27,3:25)]
View(geneInfo0)

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Sex));
geneInfo = geneInfo0[geneOrder, ]
#This data frame can be written into a text-format spreadsheet, for example by
getwd()

write.csv(geneInfo, file = "rna_WGCNA_Cell_geneInfo_sft_pwr_9_h0.2.csv")

#let's pick columns with Gene_ID, ENTREZID, and moduleColor
rna_WGCNA_Cell_Gene_ID_moduleColor_sft_pwr_9_h0.2<-
  geneInfo[,c("Gene_ID","ENTREZID","moduleColor")]
write.csv(rna_WGCNA_Cell_Gene_ID_moduleColor_sft_pwr_9_h0.2, file = "rna_WGCNA_Cell_Gene_ID_moduleColor_sft_pwr_9_h0.2.csv")

collectGarbage()

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute"))
BiocManager::install("WGCNA")
install.packages("BiocManager")
library(WGCNA)
library(Seurat)
library(flashClust)
library(cluster)
library(DOSE)
library(gplots)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#4 Interfacing network analysis with other data such as functional annotation and gene ontology
#Our previous analysis has identified several modules (labeled brown, red, and salmon) that are highly associated with weight. 
#To facilitate a biological interpretation, we would like to know the gene ontologies of the genes in the modules, 
#whether they are significantly enriched in certain functional categories etc.

#4.a Output gene lists for use with online software and services

#One option is to simply export a list of gene identifiers that can be used as input for several popular gene ontology
#and functional enrichment analysis suites such as David or AmiGO. 
#For example, we write out the LocusLinkID (entrez) codes for the brown module into a file:

View(datExprCell)
# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExprCell)
#probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = geneInfo$ENTREZID;
# $ Choose interesting modules
intModules = c("blue", "brown", "yellow", "green", "turquoise")
for (module in intModules)
{
  # Select module probes
  modGenes = (dynamicColorsCell==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
#probes = names(datExprCell)
#probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
#allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
#intModules = c("brown", "red", "salmon")
#for (module in intModules)
#{
# Select module probes
#  modGenes = (moduleColors==module)
# Get their entrez ID codes
#  modLLIDs = allLLIDs[modGenes];
# Write them into a file
#  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
#  write.table(as.data.frame(modLLIDs), file = fileName,
#              row.names = FALSE, col.names = FALSE)
#}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("Cellgenes-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

GOenr = GOenrichmentAnalysis(moduleColorsCell, allLLIDs, organism = "mouse", nBestP = 10);

tab = GOenr$bestPTerms[[4]]$enrichment

names(tab)

write.table(tab, file = "rna_WGCNA_Cell_SexTrait_GOEnrichmentTable_sft_pwr_9_h0.2.csv", sep = ",", quote = TRUE, row.names = FALSE)

keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab
View(screenTab)
View(tab)
write.table(screenTab, file = "rna_WGCNA_Cell_SexTrait_GOEnrichmentTable_shortened_sft_pwr_9_h0.2.csv", sep = ",", quote = TRUE, row.names = FALSE)


#previously:
#save(MEs, moduleLabels, moduleColors, geneTree,
#     file = "FemaleLiver-02-networkConstruction-auto.RData")
##be performed automatically:
#merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)  
#The variable merge contains various information; we will need the following:
# Numeric module labels
#moduleLabels = merge$colors;
# Convert labels to colors
#moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
#consMEs = merge$newMEs;
#
## Call an automatic merging function
#mergeEV = mergeCloseModules(datExprEV, dynamicColorsEV, cutHeight = MEDissThres, verbose = 3)
##The variable merge contains various information; we will need the following:
# Numeric module labels
#moduleLabelsEV = mergeEV$colors;
# Convert labels to colors
#moduleColorsEV = labels2colors(moduleLabelsEV)
# Eigengenes of the new merged modules:
#mergedMEsEV = mergeEV$newMEs;

table(dynamicColorsCell)
table(moduleColors)
table(dynamicLabels)
table(moduleLabels)
table(mergeCell$colors)
table(labels2colors(moduleLabelsCell))
#
save(mergedMEsCell, moduleLabelsCell, moduleColorsCell, geneTreeCell,
     file = "Cell-02-networkConstruction-man.RData")

#save(MEs, dynamicModsCell, dynamicColorsCell, geneTreeCell,
#     file = "Cell-02-networkConstruction-auto.RData")


getwd()
save.image("/scratch/user/daehyukchung/Transcriptomic_WGCNA_v4.RData")
