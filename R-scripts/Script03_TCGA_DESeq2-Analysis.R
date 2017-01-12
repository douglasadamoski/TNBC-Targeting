#
# Script03_TCGA_DESeq2-Analysis.R
#
# This R-script is a part of:
#
# Guanylate-binding protein-1 is a potential new therapeutic target for triple negative breast cancer
# (paper final citation information will be here)
# 
# Authors:
# Melissa Quintero, Douglas Adamoski, Larissa Menezes dos Reis,
# Carolline Fernanda Rodrigues Ascenção, Kaliandra de Almeida Gonçalves,
# Marília Meira Dias, Marcelo Falsarella Carazzolle, Sandra Martha Gomes Dias
#
# This script intends to perform differential expression analysis using DESeq2 and TCGA data,
# comparing Triple Negative Breast Cancer (TNBC) against Non-TNBC
#
# Inputs:
#		Non-normalized gene expression values (RSEM, RPKM, etc) and TNBC assignment
#
# Outputs:
#		Normalized gene expression values and DE gene lists
#

# Load libraries
options(scipen = 999)
library(methods)
library(EBSeq)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(edgeR)
library(devtools)
library(ggplot2)
library(ggbiplot)
library(DESeq2)
library(BiocParallel)
require(annotate)
require(org.Hs.eg.db)
library(VennDiagram)
library(plot3D)

#Create Results directory
dir.create("./03_Results", showWarnings = FALSE)
dir.create("./03_Results/TCGA_TripleVsNonTriple", showWarnings = FALSE)
dir.create("./03_Results/TCGA_TripleVsNonTriple/Heatmaps", showWarnings = FALSE)
dir.create("./03_Results/TCGA_TripleVsNonTriple/PCA_Plots", showWarnings = FALSE)
dir.create("./03_Results/TCGA_TripleVsHealthy", showWarnings = FALSE)
dir.create("./03_Results/TCGA_TripleVsHealthy/Heatmaps", showWarnings = FALSE)
dir.create("./03_Results/TCGA_TripleVsHealthy/PCA_Plots", showWarnings = FALSE)
dir.create("./03_Results/GeneralPlots", showWarnings = FALSE)
dir.create("./03_Results/GeneralPlots/Heatmaps", showWarnings = FALSE)
dir.create("./03_Results/GeneralPlots/Boxplots_AllGenes", showWarnings = FALSE)
dir.create("./03_Results/GeneralPlots/Boxplots_DEGenes", showWarnings = FALSE)
dir.create("./03_Results/GeneralPlots/PCA_Plots", showWarnings = FALSE)
dir.create("./03_Results/GeneralPlots/VennPlots", showWarnings = FALSE)

#Drop SessionInfo about packages version
sink("./03_Results/sessionInfo.txt")
sessionInfo()
print("Date executed:")
print(Sys.time())
sink()

# Reads the count table
# The list must contain genes on rows and patients on cols
# Expression values should be raw (or Raw RSEMs, unormalized)
GeneExpression_primarysite <- data.frame(read.table("./01_Results/BRCA.PrimaryTumor.txt", stringsAsFactors = FALSE))
#Change dots for dashes again
colnames(GeneExpression_primarysite) <- gsub("\\.","-", colnames(GeneExpression_primarysite))

# Reads the normal tissue
GeneExpression_normaltissue <- data.frame(read.table("./01_Results/BRCA.NormalTissue.txt", stringsAsFactors = FALSE))
#Change dots for dashes again
colnames(GeneExpression_normaltissue) <- gsub("\\.","-", colnames(GeneExpression_normaltissue))

# Define analysis groups
#Get triple vs NonTriple Information
colData <- read.table("./02_Results/DESeq2_Groups.txt", stringsAsFactors = FALSE)
data.frame(colData) -> colData

#Creates the intersection
intersect(colnames(GeneExpression_primarysite), row.names(colData)) -> conjuntocomum

#Create the table with information
GeneralInfoAll <- data.frame(matrix(ncol=1, nrow=(length(conjuntocomum)+dim(GeneExpression_normaltissue)[2]))) 
rownames(GeneralInfoAll) <- c(conjuntocomum, colnames(GeneExpression_normaltissue))
colnames(GeneralInfoAll) <- "conditions"
GeneralInfoAll$conditions <- c(colData[conjuntocomum,], rep("Healthy",dim(GeneExpression_normaltissue)[2]))

# creates a single merged table
GeneExpression_allAvailable <- cbind(GeneExpression_primarysite[,conjuntocomum], GeneExpression_normaltissue)

# It will generate an vector for normalization
NormalizacaoBullard <- QuantileNorm(GeneExpression_allAvailable,.75)

# Get the normalization values and apply to the list
GeneExpression_allAvailable_normalized <- GetNormalizedMat(GeneExpression_allAvailable, NormalizacaoBullard)

# Saves the raw complete gene expression values
write.table(GeneExpression_allAvailable,
            file="./03_Results/BRCA.all.txt")

# Saves the normalized complete gene expression values
write.table(GeneExpression_allAvailable_normalized,
            file="./03_Results/BRCA.all.normalized.txt")

# Saves the complete information
write.table(GeneralInfoAll,
            file="./03_Results/BRCA.all.information.txt")

GeneralInfoAll <- read.table(file="./03_Results/BRCA.all.information.txt", stringsAsFactors = FALSE)

######## PRIMARY SITE:
#### Triple vs Non-Triple
#### Start

#Copy the table
GeneExpression_primarysite -> GeneExpression

# Define analysis groups
colData <- read.table("./02_Results/DESeq2_Groups.txt")

#Set up it as factor and order it
colData$conditions <- as.factor(colData$conditions)
colData$conditions <- factor(colData$conditions, levels = c("Triple","NonTriple"))

# Match the Group file from Script 1 to RSEMs lists
# by order and by content
intersect(colnames(GeneExpression), row.names(colData)) -> conjuntocomum

# Create the common list
match(conjuntocomum, colnames(GeneExpression)) -> matchconjuntocomum

# Create new GeneExpression table just with matched patients
GeneExpression[,matchconjuntocomum] -> GeneExpression

# Proccedes with Bullard upper quantile normalization
# See 
#   Bullard, J. H., Purdom, E., Hansen, K. D., & Dudoit, S. (2010).
#   Evaluation of statistical methods for normalization
#   and differential expression in mRNA-Seq experiments.
#   BMC Bioinformatics, 11, 94. http://doi.org/10.1186/1471-2105-11-94
#
# It will generate an vector for normalization
NormalizacaoBullard <- QuantileNorm(GeneExpression,.75)

# Get the normalization values and apply to the list
NormalizedGeneExpression <- GetNormalizedMat(GeneExpression, NormalizacaoBullard)

# Saves normalized 
write.table(NormalizedGeneExpression, file="./03_Results/TCGA_TripleVsNonTriple/BRCA_TripleVsNonTriple_normalized.txt")

# Round values for DESeq2
# 
GeneExpression <- round(GeneExpression)

# Creates DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = GeneExpression,
                              colData = colData,
                              design = ~ conditions)

# Pre filtering for no expression genes
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Relevel
dds$conditions <- relevel(dds$conditions, ref="NonTriple")

# Perform Wald test steps
# Estimate size factors
dds <- estimateSizeFactors(dds, type="ratio", locfunc = stats::median)

# Estimate dispersions
dds <- estimateDispersions(dds, fitType = "parametric", maxit = 100, quiet = FALSE)

# Perform Wald Test
dds <- nbinomWaldTest(dds, maxit = 100, useOptim = TRUE,
                      quiet = FALSE, useT = FALSE, useQR = TRUE)

# Grab the results
resultados <- results(dds)

# Print the maplot
png(file="./03_Results/TCGA_TripleVsNonTriple/MAplot.png", res=500, width = 2000, height = 1000)
plotMA(resultados, main="DESeq2", ylim=c(-2,2))
dev.off()

# Write and read it again
write.table(resultados, file="./03_Results/TCGA_TripleVsNonTriple/TvsNT_results_DESeq2.txt")
resultados <- data.frame(read.table(file="./03_Results/TCGA_TripleVsNonTriple/TvsNT_results_DESeq2.txt", stringsAsFactors = FALSE))

# Plotting start ------
# Create table just with DE and with FC cutoff
resultadosDE <- resultados[which(resultados$padj<.05 & resultados$log2FoldChange>1),]
resultadosDE <- rbind(resultadosDE, resultados[which(resultados$padj<.05 & resultados$log2FoldChange< -1),])

# Get normalized values for heatmap
ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosDE), rownames(NormalizedGeneExpression)),]

# Create color ramp
RampaDeCor<- c(colorRampPalette(c("green","green","green","black"))(50),colorRampPalette(c("black","red","red","red"))(50))

# Create the coldata to color the conditions on heatmap
colData -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

#Select colors
ann_colors <- list(conditions = c(Triple = "#00BFC4", NonTriple = "#F8766D"))

#Heatmap generation
png(file=paste("./03_Results/TCGA_TripleVsNonTriple/Heatmaps/Correlation_Pearson_Complete_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
pheatmap(log2(ParaHeatmaps+1),
         clustering_distance_cols="correlation",
         clustering_distance_rows="correlation",
         clustering_method="complete",
         color=RampaDeCor,
         border_color=NA,
         scale="row",
         legend=TRUE,
         annotation_col=colDataHeatmap,
         annotation_colors = ann_colors)
dev.off()


#Heatmap generation
png(file=paste("./03_Results/TCGA_TripleVsNonTriple/Heatmaps/Correlation_euclidean_Complete_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
pheatmap(log2(ParaHeatmaps+1),
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method="complete",
         color=RampaDeCor,
         border_color=NA,
         scale="row",
         legend=TRUE,
         annotation_col=colDataHeatmap,
         annotation_colors = ann_colors)
dev.off()


#Heatmap generation
png(file=paste("./03_Results/TCGA_TripleVsNonTriple/Heatmaps/Correlation_Pearson_Average_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
pheatmap(log2(ParaHeatmaps+1),
         clustering_distance_cols="correlation",
         clustering_distance_rows="correlation",
         clustering_method="average",
         color=RampaDeCor,
         border_color=NA,
         scale="row",
         legend=TRUE,
         annotation_col=colDataHeatmap,
         annotation_colors = ann_colors)
dev.off()

# Start PCA Analysis
#prepara um vetor de cores e o vetor com condições
gsub("Triple", "T", colData$condition) -> ParaGrafico
gsub("NonT", "N", ParaGrafico) -> ParaGrafico
cores <- ifelse(ParaGrafico=="N","#377EB8","#E41A1C")

# Create an table just with top DE genes
GENETOP <- c("All")
resultadosDE -> resultadosPlot

# Heatmap collection
ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosPlot), rownames(NormalizedGeneExpression)),]

# Create colData
colData -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                 center = TRUE,
                 scale. = TRUE) 
ir.species <- ParaGrafico

png(file=paste("./03_Results/TCGA_TripleVsNonTriple/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
ggbiplot(ir.pca, choices = 1:2,
         obs.scale = 1,
         var.scale = 1,
         groups = ir.species, 
         ellipse = TRUE,
         circle = TRUE,
         alpha = 0.5,
         var.axes = FALSE) +
  scale_color_discrete(name = '')
#  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

# Create an table just with top DE genes
GENETOP <- 50
resultadosDE <- resultadosDE[order(resultadosDE$log2FoldChange),]
resultadosDE_Down100 <- resultadosDE[1:GENETOP,]
resultadosDE <- resultadosDE[order(-resultadosDE$log2FoldChange),]
resultadosDE_Up100 <- resultadosDE[1:GENETOP,]
rbind(resultadosDE_Down100, resultadosDE_Up100) -> resultadosPlot

ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosPlot), rownames(NormalizedGeneExpression)),]
colData -> colDataHeatmap

1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                 center = TRUE,
                 scale. = TRUE) 
ir.species <- ParaGrafico

png(file=paste("./03_Results/TCGA_TripleVsNonTriple/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
ggbiplot(ir.pca, choices = 1:2,
         obs.scale = 1,
         var.scale = 1,
         groups = ir.species, 
         ellipse = TRUE,
         circle = TRUE,
         alpha = 0.5,
         var.axes = FALSE) +
  scale_color_discrete(name = '')
#  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

# Create an table just with top DE genes
GENETOP <- 250
resultadosDE <- resultadosDE[order(resultadosDE$log2FoldChange),]
resultadosDE_Down100 <- resultadosDE[1:GENETOP,]
resultadosDE <- resultadosDE[order(-resultadosDE$log2FoldChange),]
resultadosDE_Up100 <- resultadosDE[1:GENETOP,]
rbind(resultadosDE_Down100, resultadosDE_Up100) -> resultadosPlot

ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosPlot), rownames(NormalizedGeneExpression)),]
colData -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                 center = TRUE,
                 scale. = TRUE) 
ir.species <- ParaGrafico

png(file=paste("./03_Results/TCGA_TripleVsNonTriple/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
ggbiplot(ir.pca, choices = 1:2,
         obs.scale = 1,
         var.scale = 1,
         groups = ir.species, 
         ellipse = TRUE,
         circle = TRUE,
         alpha = 0.5,
         var.axes = FALSE) +
  scale_color_discrete(name = '')
#  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

# Create an table just with top DE genes
GENETOP <- 1000
resultadosDE <- resultadosDE[order(resultadosDE$log2FoldChange),]
resultadosDE_Down100 <- resultadosDE[1:GENETOP,]
resultadosDE <- resultadosDE[order(-resultadosDE$log2FoldChange),]
resultadosDE_Up100 <- resultadosDE[1:GENETOP,]
rbind(resultadosDE_Down100, resultadosDE_Up100) -> resultadosPlot

ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosPlot), rownames(NormalizedGeneExpression)),]
colData -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                 center = TRUE,
                 scale. = TRUE) 
ir.species <- ParaGrafico

png(file=paste("./03_Results/TCGA_TripleVsNonTriple/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
ggbiplot(ir.pca, choices = 1:2,
         obs.scale = 1,
         var.scale = 1,
         groups = ir.species, 
         ellipse = TRUE,
         circle = TRUE,
         alpha = 0.5,
         var.axes = FALSE) +
  scale_color_discrete(name = '')
#  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

# START VOLCANO PLOT
resultados$Colour=rgb(100,100,100,50,maxColorValue=255)

# Set new column values to appropriate colours
#resultados$Colour[resultados$log2FoldChange > 0 & resultados$padj<0.05] <- rgb(224,162,162,50,maxColorValue=255)
resultados$Colour[resultados$log2FoldChange > 1 & resultados$padj<0.05] <- rgb(222,22,22,50,maxColorValue=255)

#resultados$Colour[resultados$log2FoldChange < 0 & resultados$padj<0.05] <- rgb(176,174,245,50,maxColorValue=255)
resultados$Colour[resultados$log2FoldChange < -1 & resultados$padj<0.05] <- rgb(56,50,237,50,maxColorValue=255)

#resultados$Colour[resultados$padj>0.95] <- rgb(25,25,25,50,maxColorValue=255)

####Volcano Plot
# Define X-axis limits
axislimits_x <- 13

# Define Y-axis limits
SmallestP <- min(resultados$padj[which(resultados$padj > 0)], na.rm=TRUE)
axislimits_y <- ceiling(-log10(SmallestP))+0.05

# Adjust pvalues to avoid NAs
VolcanoPValues <- -log10(resultados$padj + SmallestP)
VolcanoPValues[which(is.na(VolcanoPValues))] <- 0

# Simple Volcano Plot
png(file=paste("./03_Results/TCGA_TripleVsNonTriple/","VolcanoPlot_Basic_TCGA_TvsTN.png", sep=""), width=2000, height=1700, res = 300)
par(mar=c(4,4,3,2), mgp=c(2,.7,0), tck=-.01)
plot(resultados$log2FoldChange, VolcanoPValues,
     xlim=c(-axislimits_x, axislimits_x), ylim=c(0, axislimits_y),
     xlab="log2 Fold Change", ylab="-log10 p-value",
     main="Volcano Plot", cex=3.5, cex.lab=2, cex.axis=2, cex.main=1, cex.sub=2,
     pch=16, col=resultados$Colour)

abline(v = log2(2), lwd=6, col=rgb(222,22,22,150,maxColorValue=255))
abline(v = -log2(2), lwd=6, col=rgb(56,50,237,150,maxColorValue=255))

text(-12, 100, labels=
       length(which(resultados$log2FoldChange < 0 & resultados$padj < 0.05)),
     cex = 1)
text(-12, 200, labels=
       length(which(resultados$log2FoldChange < -1 & resultados$padj < 0.05)),
     cex = 1)

text(12, 100, labels=
       length(which(resultados$log2FoldChange > 0 & resultados$padj < 0.05)),
     cex = 1)
text(12, 200, labels=
       length(which(resultados$log2FoldChange > 1 & resultados$padj < 0.05)),
     cex = 1)

dev.off()

################ PRIMARY SITE vs NORMAL
#Set-up conditions and order it
Grupos <- as.factor(c(rep("Triple", length(rownames(colData)[which(colData$conditions == "Triple")])),
                      rep("Healthy", length(colnames(GeneExpression_normaltissue)))))
Grupos <- factor(Grupos, levels = c("Triple","Healthy"))

# Select Triple Patients
GeneExpression_primarysite[,rownames(colData)[which(colData$conditions == "Triple")]] -> GeneExpression

#Add NormalTissue
GeneExpression <- cbind(GeneExpression, GeneExpression_normaltissue)

# create colData object for DESeq2
colData <- matrix(ncol=1, nrow=length(Grupos))
rownames(colData) <- as.character(colnames(GeneExpression))
colnames(colData) <- "conditions"
colData[,1] <- as.character(Grupos)
colData <- data.frame(colData)

# Proccedes with Bullard upper quantile normalization for heatmap purposes
# See 
#   Bullard, J. H., Purdom, E., Hansen, K. D., & Dudoit, S. (2010).
#   Evaluation of statistical methods for normalization
#   and differential expression in mRNA-Seq experiments.
#   BMC Bioinformatics, 11, 94. http://doi.org/10.1186/1471-2105-11-94
#
# It will generate an vector for normalization
NormalizacaoBullard <- QuantileNorm(GeneExpression,.75)

# Get the normalization values and apply to the list
NormalizedGeneExpression <- GetNormalizedMat(GeneExpression, NormalizacaoBullard)

# Saves normalized 
write.table(NormalizedGeneExpression, file="./03_Results/TCGA_TripleVsHealthy/BRCA_TripleVsHealthy_normalized.txt")


# Round values for DESeq2
# 
GeneExpression <- round(GeneExpression)

# Creates DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = GeneExpression,
                              colData = colData,
                              design = ~ conditions)

# Pre filtering for no expression genes
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Relevel
dds$conditions <- relevel(dds$conditions, ref="Normal")

# Perform Wald test steps
# Estimate size factors
dds <- estimateSizeFactors(dds, type="ratio", locfunc = stats::median)

# Estimate dispersions
dds <- estimateDispersions(dds, fitType = "parametric", maxit = 100, quiet = FALSE)

# Perform Wald Test
dds <- nbinomWaldTest(dds, maxit = 100, useOptim = TRUE,
                      quiet = FALSE, useT = FALSE, useQR = TRUE)

# Grab the results
resultados <- results(dds)

# Print the maplot
png(file="./03_Results/TCGA_TripleVsHealthy/MAplot.png", res=500, width = 2000, height = 1000)
plotMA(resultados, main="DESeq2", ylim=c(-2,2))
dev.off()

# Write it
write.table(resultados, file="./03_Results/TCGA_TripleVsHealthy/TvsH_results_DESeq2.txt")
resultados <- data.frame(read.table(file="./03_Results/TCGA_TripleVsHealthy/TvsH_results_DESeq2.txt", stringsAsFactors = FALSE))

# Plotting start ------
# Create table just with DE and with FC cutoff
resultadosDE <- resultados[which(resultados$padj<.05 & resultados$log2FoldChange>1),]
resultadosDE <- rbind(resultadosDE, resultados[which(resultados$padj<.05 & resultados$log2FoldChange< -1),])

# Get normalized values for heatmap
ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosDE), rownames(NormalizedGeneExpression)),]

# Create color ramp
RampaDeCor<- c(colorRampPalette(c("green","green","green","black"))(50),colorRampPalette(c("black","red","red","red"))(50))

# Create the coldata to color the conditions on heatmap
colData -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

# Select colors
ann_colors <- list(conditions = c(Triple = "#00BFC4", Healthy = "#ec7014"))

# Heatmap generation
png(file=paste("./03_Results/TCGA_TripleVsHealthy/Heatmaps/Correlation_Pearson_Complete_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
pheatmap(log2(ParaHeatmaps+1),
         clustering_distance_cols="correlation",
         clustering_distance_rows="correlation",
         clustering_method="complete",
         color=RampaDeCor,
         border_color=NA,
         scale="row",
         legend=TRUE,
         annotation_col=colDataHeatmap,
         annotation_colors = ann_colors)
dev.off()


# Heatmap generation
png(file=paste("./03_Results/TCGA_TripleVsHealthy/Heatmaps/Correlation_euclidean_Complete_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
pheatmap(log2(ParaHeatmaps+1),
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method="complete",
         color=RampaDeCor,
         border_color=NA,
         scale="row",
         legend=TRUE,
         annotation_col=colDataHeatmap,
         annotation_colors = ann_colors)
dev.off()


# Heatmap generation
png(file=paste("./03_Results/TCGA_TripleVsHealthy/Heatmaps/Correlation_Pearson_Average_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
pheatmap(log2(ParaHeatmaps+1),
         clustering_distance_cols="correlation",
         clustering_distance_rows="correlation",
         clustering_method="average",
         color=RampaDeCor,
         border_color=NA,
         scale="row",
         legend=TRUE,
         annotation_col=colDataHeatmap,
         annotation_colors = ann_colors)
dev.off()

# Start PCA Analysis
# create the condition vector
gsub("Triple", "T", colDataHeatmap$conditions) -> ParaGrafico
gsub("Normal", "H", ParaGrafico) -> ParaGrafico
cores <- ifelse(ParaGrafico=="H","#ec7014","#E41A1C")

# DE Table
GENETOP <- c("All")
resultadosDE -> resultadosPlot

ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosPlot), rownames(NormalizedGeneExpression)),]
colData -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                 center = TRUE,
                 scale. = TRUE) 
ir.species <- ParaGrafico

png(file=paste("./03_Results/TCGA_TripleVsHealthy/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
ggbiplot(ir.pca, choices = 1:2,
         obs.scale = 1,
         var.scale = 1,
         groups = ir.species, 
         ellipse = TRUE,
         circle = TRUE,
         alpha = 0.5,
         var.axes = FALSE) +
  scale_color_discrete(name = '')
#  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

# DE Table
GENETOP <- 50
resultadosDE <- resultadosDE[order(resultadosDE$log2FoldChange),]
resultadosDE_Down100 <- resultadosDE[1:GENETOP,]
resultadosDE <- resultadosDE[order(-resultadosDE$log2FoldChange),]
resultadosDE_Up100 <- resultadosDE[1:GENETOP,]
rbind(resultadosDE_Down100, resultadosDE_Up100) -> resultadosPlot


ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosPlot), rownames(NormalizedGeneExpression)),]
colData -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                 center = TRUE,
                 scale. = TRUE) 
ir.species <- ParaGrafico

png(file=paste("./03_Results/TCGA_TripleVsHealthy/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
ggbiplot(ir.pca, choices = 1:2,
         obs.scale = 1,
         var.scale = 1,
         groups = ir.species, 
         ellipse = TRUE,
         circle = TRUE,
         alpha = 0.5,
         var.axes = FALSE) +
  scale_color_discrete(name = '')
#  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

# DE Table
GENETOP <- 250
resultadosDE <- resultadosDE[order(resultadosDE$log2FoldChange),]
resultadosDE_Down100 <- resultadosDE[1:GENETOP,]
resultadosDE <- resultadosDE[order(-resultadosDE$log2FoldChange),]
resultadosDE_Up100 <- resultadosDE[1:GENETOP,]
rbind(resultadosDE_Down100, resultadosDE_Up100) -> resultadosPlot

ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosPlot), rownames(NormalizedGeneExpression)),]
colData -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                 center = TRUE,
                 scale. = TRUE) 
ir.species <- ParaGrafico

png(file=paste("./03_Results/TCGA_TripleVsHealthy/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
ggbiplot(ir.pca, choices = 1:2,
         obs.scale = 1,
         var.scale = 1,
         groups = ir.species, 
         ellipse = TRUE,
         circle = TRUE,
         alpha = 0.5,
         var.axes = FALSE) +
  scale_color_discrete(name = '')
#  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

# DE Table
GENETOP <- 1000
resultadosDE <- resultadosDE[order(resultadosDE$log2FoldChange),]
resultadosDE_Down100 <- resultadosDE[1:GENETOP,]
resultadosDE <- resultadosDE[order(-resultadosDE$log2FoldChange),]
resultadosDE_Up100 <- resultadosDE[1:GENETOP,]
rbind(resultadosDE_Down100, resultadosDE_Up100) -> resultadosPlot

ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosPlot), rownames(NormalizedGeneExpression)),]
colData -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                 center = TRUE,
                 scale. = TRUE) 
ir.species <- ParaGrafico

png(file=paste("./03_Results/TCGA_TripleVsHealthy/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
ggbiplot(ir.pca, choices = 1:2,
         obs.scale = 1,
         var.scale = 1,
         groups = ir.species, 
         ellipse = TRUE,
         circle = TRUE,
         alpha = 0.5,
         var.axes = FALSE) +
  scale_color_discrete(name = '')
#  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

# START VOLCANO PLOT
resultados$Colour=rgb(100,100,100,50,maxColorValue=255)

# Set new column values to appropriate colours
resultados$Colour[resultados$log2FoldChange > 1 & resultados$padj<0.05] <- rgb(222,22,22,50,maxColorValue=255)
resultados$Colour[resultados$log2FoldChange < -1 & resultados$padj<0.05] <- rgb(56,50,237,50,maxColorValue=255)

####Volcano Plot
# Define X-axis limits
axislimits_x <- 13

# Define Y-axis limits
SmallestP <- min(resultados$padj[which(resultados$padj > 0)], na.rm=TRUE)
axislimits_y <- ceiling(-log10(SmallestP))+0.05

# Adjust pvalues to avoid NAs
VolcanoPValues <- -log10(resultados$padj + SmallestP)
VolcanoPValues[which(is.na(VolcanoPValues))] <- 0

# Simple Volcano Plot
png(file=paste("./03_Results/TCGA_TripleVsHealthy/","VolcanoPlot_Basic_TCGA_TvsTN.png", sep=""), width=2000, height=1700, res = 300)
par(mar=c(4,4,3,2), mgp=c(2,.7,0), tck=-.01)
plot(resultados$log2FoldChange, VolcanoPValues,
     xlim=c(-axislimits_x, axislimits_x), ylim=c(0, axislimits_y),
     xlab="log2 Fold Change", ylab="-log10 p-value",
     main="Volcano Plot", cex=3.5, cex.lab=2, cex.axis=2, cex.main=1, cex.sub=2,
     pch=16, col=resultados$Colour)

abline(v = log2(2), lwd=6, col=rgb(222,22,22,150,maxColorValue=255))
abline(v = -log2(2), lwd=6, col=rgb(56,50,237,150,maxColorValue=255))

text(-12, 100, labels=
       length(which(resultados$log2FoldChange < 0 & resultados$padj < 0.05)),
     cex = 1)
text(-12, 200, labels=
       length(which(resultados$log2FoldChange < -1 & resultados$padj < 0.05)),
     cex = 1)

text(12, 100, labels=
       length(which(resultados$log2FoldChange > 0 & resultados$padj < 0.05)),
     cex = 1)
text(12, 200, labels=
       length(which(resultados$log2FoldChange > 1 & resultados$padj < 0.05)),
     cex = 1)

dev.off()

### START INTERSECTION BETWEEN T vs NT and TN vs HEALTHY
# Read results back
resultadosTvsNT <- data.frame(read.table(file="./03_Results/TCGA_TripleVsNonTriple/TvsNT_results_DESeq2.txt", stringsAsFactors = FALSE))
resultadosTvsH <- data.frame(read.table(file="./03_Results/TCGA_TripleVsHealthy/TvsH_results_DESeq2.txt", stringsAsFactors = FALSE))

#Change colnames in order to cbind
colnames(resultadosTvsNT) <- paste(colnames(resultadosTvsNT),"_TvsNT", sep="")
colnames(resultadosTvsH) <- paste(colnames(resultadosTvsH),"_TvsH", sep="")

#Intersect the list
intersection <- intersect(rownames(resultadosTvsNT),
rownames(resultadosTvsH))

# Create the intersected result table
ResultsIntersected <- cbind(resultadosTvsNT[intersection,],
resultadosTvsH[intersection,])

#Create the both significant column
ResultsIntersected$BothTarget <- "No"

# Defines if the target is DE vsNonTriple, vsHealthy or both
ResultsIntersected$BothTarget[
  which(ResultsIntersected$padj_TvsNT < 0.05 & abs(ResultsIntersected$log2FoldChange_TvsNT) > 1)
] <- "vsNonTriple"

ResultsIntersected$BothTarget[
  which(ResultsIntersected$padj_TvsH < 0.05 & abs(ResultsIntersected$log2FoldChange_TvsH) > 1)
] <- "vsHealthy"       
        
ResultsIntersected$BothTarget[
  which(
    ResultsIntersected$padj_TvsNT < 0.05 & ResultsIntersected$log2FoldChange_TvsNT > 1 &
    ResultsIntersected$padj_TvsH < 0.05 & ResultsIntersected$log2FoldChange_TvsH > 1)
  ] <- "Both" 

ResultsIntersected$BothTarget[
  which(
    ResultsIntersected$padj_TvsNT < 0.05 & ResultsIntersected$log2FoldChange_TvsNT < -1 &
      ResultsIntersected$padj_TvsH < 0.05 & ResultsIntersected$log2FoldChange_TvsH < -1)
  ] <- "Both" 

# Write and read it again
write.table(ResultsIntersected, file="./03_Results/results_DESeq2_intersected.txt")
write.csv2(ResultsIntersected, file="./03_Results/results_DESeq2_intersected.csv")
ResultsIntersected <- data.frame(read.table(file="./03_Results/results_DESeq2_intersected.txt", stringsAsFactors = FALSE))


# Plotting start ------
### VennDiagramns
# Defines if the target is DE vsNonTriple, vsHealthy or both

TvsNT_up <- rownames(ResultsIntersected)[
which(ResultsIntersected$padj_TvsNT < 0.05 & ResultsIntersected$log2FoldChange_TvsNT > 1)]

TvsNT_down <- rownames(ResultsIntersected)[
  which(ResultsIntersected$padj_TvsNT < 0.05 & ResultsIntersected$log2FoldChange_TvsNT < -1)]

TvsH_up <- rownames(ResultsIntersected)[
which(ResultsIntersected$padj_TvsH < 0.05 & ResultsIntersected$log2FoldChange_TvsH > 1)]

TvsH_down <- rownames(ResultsIntersected)[
  which(ResultsIntersected$padj_TvsH < 0.05 & ResultsIntersected$log2FoldChange_TvsH < -1)]

# All DE List
venn.plot <- venn.diagram(
  list(A = c(TvsNT_up,TvsNT_down), B = c(TvsH_up,TvsH_down)),
  height = 3000, width = 3000, resolution = 500, imagetype = "png",
  main="All DE List TCGA",
  main.cex=1,
  print.mode="raw",
  lwd=c("4","4"),
  col=brewer.pal(8,'Set3')[c(7,6)],
  fill=brewer.pal(8,'Set3')[c(7,6)],
  alpha=0.5,
  cat.cex = 1.5,
  cex=c(2,2,2),
  margin = 0.12,
  cat.dist=c(0.03,0.04),
  category.names = c("TvsNT","TvsH"),
  filename="./03_Results/GeneralPlots/VennPlots/VennAll_Both.png")

# UP Regulated
venn.plot <- venn.diagram(
  list(A = c(TvsNT_up), B = c(TvsH_up)),
  height = 3000, width = 3000, resolution = 500, imagetype = "png",
  main="UP List TCGA",
  main.cex=1,
  print.mode="raw",
  lwd=c("4","4"),
  col=c("#F8766D"),
  fill=c("#F8766D"),
  alpha=0.5,
  cat.cex = 1.5,
  cex=c(2,2,2),
  margin = 0.12,
  cat.dist=c(0.03,0.04),
  category.names = c("TvsNT","TvsH"),
  filename="./03_Results/GeneralPlots/VennPlots/VennAll_Up.png")

# Down Regulated
venn.plot <- venn.diagram(
  list(A = c(TvsNT_down), B = c(TvsH_down)),
  height = 3000, width = 3000, resolution = 500, imagetype = "png",
  main="Down List TCGA",
  main.cex=1,
  print.mode="raw",
  lwd=c("4","4"),
  col=c("#00BFC4"),
  fill=c("#00BFC4"),
  alpha=0.5,
  cat.cex = 1.5,
  cex=c(2,2,2),
  margin = 0.12,
  cat.dist=c(0.03,0.04),
  category.names = c("TvsNT","TvsH"),
  filename="./03_Results/GeneralPlots/VennPlots/VennAll_Down.png")

#Rename the normalized matrix
NormalizedGeneExpression <- GeneExpression_allAvailable_normalized

# Create table just with DE and with FC cutoff
resultadosDE <- ResultsIntersected[which(ResultsIntersected$BothTarget == "Both"),]

# Get normalized values for heatmap
ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosDE), rownames(NormalizedGeneExpression)),]

# Create color ramp
RampaDeCor<- c(colorRampPalette(c("green","green","green","black"))(50),colorRampPalette(c("black","red","red","red"))(50))

# Create the coldata to color the conditions on heatmap
GeneralInfoAll -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)

#Select colors
ann_colors <- list(conditions = c(Triple = "#00BFC4",  NonTriple = "#3AC664", Healthy = "#F8766D"))

#Heatmap generation
png(file=paste("./03_Results/GeneralPlots/Heatmaps/Correlation_Pearson_Complete_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
pheatmap(log2(ParaHeatmaps+1),
         clustering_distance_cols="correlation",
         clustering_distance_rows="correlation",
         clustering_method="complete",
         color=RampaDeCor,
         border_color=NA,
         scale="row",
         legend=TRUE,
         annotation_col=colDataHeatmap,
         annotation_colors = ann_colors)
dev.off()

#Heatmap generation
png(file=paste("./03_Results/GeneralPlots/Heatmaps/Correlation_euclidean_Complete_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
pheatmap(log2(ParaHeatmaps+1),
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method="complete",
         color=RampaDeCor,
         border_color=NA,
         scale="row",
         legend=TRUE,
         annotation_col=colDataHeatmap,
         annotation_colors = ann_colors)
dev.off()

#Heatmap generation
png(file=paste("./03_Results/GeneralPlots/Heatmaps/Correlation_Pearson_Average_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
pheatmap(log2(ParaHeatmaps+1),
         clustering_distance_cols="correlation",
         clustering_distance_rows="correlation",
         clustering_method="average",
         color=RampaDeCor,
         border_color=NA,
         scale="row",
         legend=TRUE,
         annotation_col=colDataHeatmap,
         annotation_colors = ann_colors)
dev.off()

# Start PCA Analysis
# create color vector
gsub("Triple", "T", colDataHeatmap$conditions) -> ParaGrafico
gsub("NonT", "N", ParaGrafico) -> ParaGrafico
gsub("Healthy", "H", ParaGrafico) -> ParaGrafico
#cores <- ifelse(ParaGrafico=="H","#ec7014","#E41A1C")

# Select all genes for analysis
GENETOP <- c("All")
resultadosDE -> resultadosPlot

# Get the normalized values for PCA plot
ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosPlot), rownames(NormalizedGeneExpression)),]

# create coldata
GeneralInfoAll -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)

ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                 center = TRUE,
                 scale. = TRUE) 
ir.species <- ParaGrafico

png(file=paste("./03_Results/GeneralPlots/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
ggbiplot(ir.pca, choices = 1:2,
         obs.scale = 1,
         var.scale = 1,
         groups = ir.species, 
         ellipse = TRUE,
         circle = TRUE,
         alpha = 0.5,
         var.axes = FALSE) +
  scale_color_discrete(name = '')
#  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

#################### Naming correction TCGA data
#################### START

# this block is necessary due naming clash between TCGAs modified hg19 against original hg19 naming annotation

#Split in two vectors
FirstName <- sapply(rownames(ResultsIntersected),
function(w) unlist(strsplit(w,"\\|"))[1])
names(FirstName) <- NULL

SecondName <- sapply(rownames(ResultsIntersected),
                    function(w) unlist(strsplit(w,"\\|"))[2])
names(SecondName) <- NULL

# Create the recipient table
TCGA_naming <- matrix(ncol=7, nrow=length(rownames(ResultsIntersected)))
rownames(TCGA_naming) <- rownames(ResultsIntersected)
colnames(TCGA_naming) <- c("UCSC_symbol_split","UCSC_ID_split","geneSymbol_fromID", "geneID_fromAlias", "geneSymbol_fromAliasID","geneSymbol_final","geneID_final")
TCGA_naming <- data.frame(TCGA_naming)

# Fill with original values
TCGA_naming$UCSC_symbol_split <- FirstName
TCGA_naming$UCSC_ID_split <- SecondName

# Get the GeneIdList
GeneIDList <- data.frame(toTable(org.Hs.egSYMBOL))
GeneAliasList <- data.frame(toTable(org.Hs.egALIAS2EG))

#
k <- 0
for(TCGANow in rownames(TCGA_naming)){
  k <- k+1
  # Use ID to find Symbol
  Position <- match(TCGA_naming$UCSC_ID_split[k], GeneIDList$gene_id)
  GeneIDList$symbol[Position] -> TCGA_naming$geneSymbol_fromID[k]
 
  # Use Alias to find ID
  Position <- match(TCGA_naming$UCSC_symbol_split[k], GeneAliasList$alias_symbol)
  GeneAliasList$gene_id[Position] -> TCGA_naming$geneID_fromAlias[k] 
  
  # Use ID from Alias to find Symbol
  Position <- match(TCGA_naming$geneID_fromAlias[k], GeneIDList$gene_id)
  GeneIDList$symbol[Position] -> TCGA_naming$geneSymbol_fromAliasID[k]

  #
  if(is.na(TCGA_naming$geneSymbol_fromID[k])){
        if(is.na(TCGA_naming$geneSymbol_fromAliasID[k])){
      TCGA_naming$UCSC_symbol_split[k] ->  TCGA_naming$geneSymbol_final[k]
      TCGA_naming$UCSC_ID_split[k] ->  TCGA_naming$geneID_final[k]
    } else {
      TCGA_naming$geneID_fromAlias[k] ->  TCGA_naming$geneID_final[k]
      TCGA_naming$geneSymbol_fromAliasID[k] ->  TCGA_naming$geneSymbol_final[k]
  
    }
  } else {
    TCGA_naming$UCSC_ID_split[k] ->  TCGA_naming$geneID_final[k]
    TCGA_naming$geneSymbol_fromID[k] ->  TCGA_naming$geneSymbol_final[k]
      }
  }

# Solve the interrogations
for(creamcheese in 1:length(TCGA_naming$geneSymbol_final)){
  if(TCGA_naming$geneSymbol_final[creamcheese]=="?"){
    TCGA_naming$geneID_final[creamcheese] -> TCGA_naming$geneSymbol_final[creamcheese]
  } else { }}

# Puts old names on duplicates
DuplicatedTCGA <- TCGA_naming$geneSymbol_final[which(duplicated(TCGA_naming$geneSymbol_final))]
#Solve the duplicates by the 'hard' way
for(carrot in 1:length(DuplicatedTCGA)){
  match(TCGA_naming$geneSymbol_final, DuplicatedTCGA[carrot]) -> TCGANowMatch
  which(!is.na(TCGANowMatch)) -> NewPositions
  for(strawberry in 1:length(NewPositions)){
    TCGA_naming$UCSC_symbol_split[NewPositions[strawberry]] -> TCGA_naming$geneSymbol_final[NewPositions[strawberry]]
    TCGA_naming$UCSC_ID_split[NewPositions[strawberry]] -> TCGA_naming$geneID_final[NewPositions[strawberry]]
  }
}

# Perform a colbind with results table
ResultsIntersected <- cbind(TCGA_naming[,c("geneSymbol_final", "geneID_final")], ResultsIntersected)

# Write and read back
write.table(TCGA_naming, file="./03_Results/Naming_intersected.txt")
write.table(ResultsIntersected, file="./03_Results/results_DESeq2_intersected_naming.txt")
write.csv2(ResultsIntersected, file="./03_Results/results_DESeq2_intersected_naming.csv")
ResultsIntersected <- data.frame(read.table(file="./03_Results/results_DESeq2_intersected_naming.txt", stringsAsFactors = FALSE))

# Create the 3D volcano plot for both comparisons
#Generate the vector with the highest pValue for a single comparison
Maximum_padj <- rep(1, dim(ResultsIntersected)[1])
for(k in 1:dim(ResultsIntersected)[1]){
  Maximum_padj[k]  <- max(ResultsIntersected$padj_TvsNT[k],ResultsIntersected$padj_TvsH[k]) 
}
#Replace NA for 1
Maximum_padj[which(is.na(Maximum_padj))] <- 1

# Creates the 3D scatterplot
# x, y and z coordinates
x <- ResultsIntersected$log2FoldChange_TvsNT
y <- ResultsIntersected$log2FoldChange_TvsH
z <- -log10(Maximum_padj)

# Create the color vector
ColorVector <- rep(0, dim(ResultsIntersected)[1])
ColorVector[which(ResultsIntersected$log2FoldChange_TvsNT > 1 &
  ResultsIntersected$log2FoldChange_TvsH > 1 &
  ResultsIntersected$BothTarget == "Both")] <- 1
ColorVector[which(ResultsIntersected$log2FoldChange_TvsNT < -1 &
                    ResultsIntersected$log2FoldChange_TvsH < -1 &
                    ResultsIntersected$BothTarget == "Both")] <- -1

# Scatterplot 3d
png(paste("./03_Results/GeneralPlots/Volcano3d",".png", sep=""), width=5000, height=7000, res=1200, type='cairo')
scatter3D(x, y, z, theta=57, phi = 20, 
          bty = "g",
          colvar=ColorVector,
          colkey = FALSE,
          col=c(rgb(56,50,237,100,maxColorValue=255), rgb(100,100,100,25,maxColorValue=255), rgb(222,22,22,100,maxColorValue=255)),
          pch = 1, cex = 0.8, ticktype = "detailed",
          xlab="NT",
          ylab="H",
          zlab="p")
dev.off()

# Boxplot genes in both lists
# This should take long
# Set the levels
mylevels <- c("Healthy", "NonTriple", "Triple")

##################################### >>>>> ONLY DE GENES
# Create table just with DE and with FC cutoff
resultadosDE <- ResultsIntersected[which(ResultsIntersected$BothTarget == "Both"),]

#Define targets
ActualGeneList <- rownames(resultadosDE)

#Empties k value
k <- 0
# Boxplot Loop START >>>>>>>>>>>>>>>>
for(Target in ActualGeneList){
  k <- k+1
  
  # Get the gene symbol
  alvos <- resultadosDE$geneSymbol_final[k]
  
  # Stores the value in our temporary table
  PlottingNow <- data.frame(matrix(nrow=dim(GeneExpression_allAvailable_normalized)[2], ncol=0))
  rownames(PlottingNow) <- colnames(GeneExpression_allAvailable_normalized)
  PlottingNow$value <- log2(GeneExpression_allAvailable_normalized[Target,]+1)
  PlottingNow$condition <- GeneralInfoAll$conditions
  
  # Perform the boxplot itself
  png(paste("./03_Results/GeneralPlots/Boxplots_DEGenes/TCGA_",alvos,".png", sep=""), width=5000, height=7000, res=1200, type='cairo')
  boxplot(PlottingNow$value~PlottingNow$condition,
          las = 1,
          col = c("#F5B4AF","#80D99A","#79D8DB"),
          border = c("#F8766D","#3AC664","#00BFC4"),
          outline=TRUE,
          range=0.5,
          boxwex=0.8,
          notch = TRUE,
          lwd=2,
          #width=levelProportions,
          ylab ="log2(RSEM+1)",
          main=alvos,
          par(cex.lab=1.5, cex.axis=1),
          outpch=NA
  )
  # Add points by jittering over boxplot
  for(i in 1:length(mylevels)){
    thislevel<-mylevels[i]
    thisvalues<-PlottingNow$value[PlottingNow$condition==thislevel]
    # take the x-axis indices and add a jitter, proportional to the N in each level
    myjitter <- jitter(rep(i, length(thisvalues)),
                     amount=0.37
                     #amount=levelProportions[i]/2
    )
    points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.1), cex=1) 
  }
  dev.off()
}
################ >>>>> ONLY DE GENES

##################################### >>>>> All GENES
# Create table just with DE and with FC cutoff
resultadosDE <- ResultsIntersected

#Define targets
ActualGeneList <- rownames(resultadosDE)

#Empties k value
k <- 0
# Boxplot Loop START >>>>>>>>>>>>>>>>
for(Target in ActualGeneList){
  k <- k+1
  
  # Get the gene symbol
  alvos <- resultadosDE$geneSymbol_final[k]
  
  # Stores the value in our temporary table
  PlottingNow <- data.frame(matrix(nrow=dim(GeneExpression_allAvailable_normalized)[2], ncol=0))
  rownames(PlottingNow) <- colnames(GeneExpression_allAvailable_normalized)
  PlottingNow$value <- as.numeric(log2(GeneExpression_allAvailable_normalized[Target,]+1))
  PlottingNow$condition <- GeneralInfoAll$conditions
  
  # perform the boxplot itself
  png(paste("./03_Results/GeneralPlots/Boxplots_AllGenes/TCGA_",alvos,".png", sep=""), width=5000, height=7000, res=1200, type='cairo')
  boxplot(PlottingNow$value~PlottingNow$condition,
          las = 1,
          col = c("#F5B4AF","#80D99A","#79D8DB"),
          border = c("#F8766D","#3AC664","#00BFC4"),
          outline=TRUE,
          range=0.5,
          boxwex=0.8,
          notch = TRUE,
          lwd=2,
          #width=levelProportions,
          ylab ="log2(RSEM+1)",
          main=alvos,
          par(cex.lab=1.5, cex.axis=1),
          outpch=NA
  )
  # Add points by jittering over boxplot
  for(i in 1:length(mylevels)){
    thislevel<-mylevels[i]
    thisvalues<-PlottingNow$value[PlottingNow$condition==thislevel]
    # take the x-axis indices and add a jitter, proportional to the N in each level
    myjitter <- jitter(rep(i, length(thisvalues)),
                       amount=0.37
                       #amount=levelProportions[i]/2
    )
    points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.1), cex=1) 
  }
  dev.off()
}
################ >>>>> ALL GENES

# Finishes the script
