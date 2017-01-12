#
# Script05_Cell_Line_DESeq2_PathwayEnrichment.R
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
# This script intends to perform DE genes analysis and pathway enrichments using RNA-seqs from cell lines
#
# Inputs:
#		Non-normalized gene expression values (RSEM, RPKM, etc) and Grouping scheme
#
# Outputs:
#		Normalized gene expression values and classification table
#		DE genes list
#		Pathway enrichments
#		A lot of cool graphs
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
library(plot3D)
require(GO.db)
library(goseq)
library(clusterProfiler)
library(DOSE)
library(ReactomePA)
library(DO.db)
library(reactome.db)
library(VennDiagram)
library(UniProt.ws)
options(bitmapType='cairo')

# Create Results directory
dir.create("./05_Results", showWarnings = FALSE)
dir.create("./05_Results/Intersection", showWarnings = FALSE)
dir.create("./05_Results/Intersection/VennDiagrams", showWarnings = FALSE)
dir.create("./05_Results/Intersection/Correlations", showWarnings = FALSE)
dir.create("./05_Results/Intersection/Pathway", showWarnings = FALSE)
dir.create("./05_Results/CellLines", showWarnings = FALSE)
dir.create("./05_Results/CellLines/Heatmaps", showWarnings = FALSE)
dir.create("./05_Results/CellLines/PCA_Plots", showWarnings = FALSE)
dir.create("./05_Results/CellLines/Boxplots_DEGenes", showWarnings = FALSE)
dir.create("./05_Results/CellLines/Boxplots_AllGenes", showWarnings = FALSE)
dir.create("./05_Results/CellLines/qPCR", showWarnings = FALSE)
dir.create("./05_Results/CellLines/qPCR_Spearman", showWarnings = FALSE)
dir.create("./05_Results/CellLines/Densities", showWarnings = FALSE)
dir.create("./05_Results/CellLines/Pathway", showWarnings = FALSE)
dir.create("./05_Results/CellLines/Pathway/GO", showWarnings = FALSE)
dir.create("./05_Results/CellLines/Pathway/REACTOME", showWarnings = FALSE)
dir.create("./05_Results/CellLines/Pathway/DOSE", showWarnings = FALSE)

#Drop SessionInfo about packages version
sink("./05_Results/sessionInfo.txt")
sessionInfo()
print("Date executed:")
print(Sys.time())
sink()

# Read classification of cells
read.table(file="./ExternalFiles/Classification_CellLines.txt") -> PrepreColData

# Generates the matrix for final colData
PrecolData <- data.frame(matrix(nrow=length(PrepreColData[,1]),ncol=2))
colnames(PrecolData) <- c("conditions","type")
rownames(PrecolData) <- rownames(PrepreColData)
factor(PrepreColData[,1], levels = c("NonTriple","Triple")) -> PrecolData[,1]
c("paired-end") -> PrecolData[,2]
PrecolData -> colData

# Creates the expression table
# Open the first file to use as model
Previa <- read.table(paste("./ExternalFiles/RSEM_CellLines/",rownames(colData)[1],".genes.results",sep=""), header = TRUE)

# Create the final table with correct number of columns
# from colData file and gene number
countData <- matrix(nrow=length(rownames(Previa)), ncol=length(rownames(colData)))

# Add names
rownames(countData) <- Previa$gene_id
colnames(countData) <- rownames(colData)

# Clear for memory purposes
remove(Previa)

# For loop to read and get everything done
i <- 0
for(kiwi in rownames(colData)){
  i + 1 -> i
  AmostraNow <- read.table(paste("./ExternalFiles/RSEM_CellLines/",kiwi,".genes.results", sep=""), header = TRUE)
  geneNamesNow <- as.character(AmostraNow$gene_id)
  rsemNow <- as.numeric(AmostraNow$expected_count)
  
  if(isTRUE(unique(geneNamesNow==rownames(countData)))){
  # Add RSEM value
    rsemNow ->countData[,i] 
  } else {
      stop("Gene lists are distinct!!!")
    }
}

# Just a name change
GeneExpression <- countData

# It will generate an vector for normalization
NormalizacaoBullard <- QuantileNorm(GeneExpression,.75)

# Get the normalization values and apply to the list
NormalizedGeneExpression <- GetNormalizedMat(GeneExpression, NormalizacaoBullard)

# Saves the raw complete gene expression values
write.table(GeneExpression,
            file="./05_Results/CellLines/RSEM_CellLines.txt")
write.table(NormalizedGeneExpression,
            file="./05_Results/CellLines/RSEM_CellLines.normalized.txt")
			
#Read it back
NormalizedGeneExpression <- read.table(file="./05_Results/CellLines/RSEM_CellLines.normalized.txt", stringsAsFactors=FALSE)

#Set up it as factor and order it
colData$conditions <- as.factor(colData$conditions)
colData$conditions <- factor(colData$conditions, levels = c("Triple","NonTriple"))

# Round values for DESeq2
# From Michael Love (DESeq2 developer) bioconductor support page (https://support.bioconductor.org/p/51577/)
# also contradicted in https://support.bioconductor.org/p/57647/
# but we decided to round-and-go
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
dds <- estimateDispersions(dds, fitType = "parametric", maxit = 200, quiet = FALSE)

# Perform Wald Test
dds <- nbinomWaldTest(dds, maxit = 100, useOptim = TRUE,
                      quiet = FALSE, useT = FALSE, useQR = TRUE)

# Grab the results
resultados <- results(dds)

# Print the maplot
png(file="./05_Results/CellLines/MAplot.png", res=500, width = 2000, height = 1000)
plotMA(resultados, main="DESeq2", ylim=c(-2,2))
dev.off()

# Write and read it again
write.table(resultados, file="./05_Results/CellLines/CellLines_TvsNT_results_DESeq2.txt")
resultados <- data.frame(read.table(file="./05_Results/CellLines/CellLines_TvsNT_results_DESeq2.txt", stringsAsFactors = FALSE))

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
#1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
#1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

#Select colors
ann_colors <- list(conditions = c(Triple = "#00BFC4", NonTriple = "#3AC664"))

#Heatmap generation
png(file=paste("./05_Results/CellLines/Heatmaps/Correlation_Pearson_Complete_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
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
png(file=paste("./05_Results/CellLines/Heatmaps/Correlation_euclidean_Complete_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
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
png(file=paste("./05_Results/CellLines/Heatmaps/Correlation_Pearson_Average_Normalized_Top_0.05.png", sep=""), res=500, width = 5000, height = 6000)
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
# Prepare condition vector
gsub("Triple", "T", colData$condition) -> ParaGrafico
gsub("NonT", "N", ParaGrafico) -> ParaGrafico

# Create the table for DE
GENETOP <- c("All")
resultadosDE -> resultadosPlot

# Heatmap separation
ParaHeatmaps <- NormalizedGeneExpression[match(rownames(resultadosPlot), rownames(NormalizedGeneExpression)),]

# Creates the colData
colData -> colDataHeatmap
1:length(rownames(colDataHeatmap)) -> rownames(colDataHeatmap)
1:length(colnames(ParaHeatmaps)) -> colnames(ParaHeatmaps)
colDataHeatmap$type <- NULL

ir.pca <- prcomp(t(log2(ParaHeatmaps+1)),
                 center = TRUE,
                 scale. = TRUE) 
ir.species <- ParaGrafico

png(file=paste("./05_Results/CellLines/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
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

# Create the table for DE
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

png(file=paste("./05_Results/CellLines/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
ggbiplot(ir.pca, choices = 1:2,
         obs.scale = 1,
         var.scale = 1,
         groups = ir.species, 
         ellipse = TRUE,
         circle = TRUE,
         alpha = 0.5,
         var.axes = FALSE) +
  scale_color_discrete(name = '')
# theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

# Create the table for DE
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

png(file=paste("./05_Results/CellLines/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
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

# Create the table for DE
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

png(file=paste("./05_Results/CellLines/PCA_Plots/Top_",GENETOP,"x2_0.05.png", sep=""), res=800, height = 5000, width=3000)
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
# Define os limites do eixo x
axislimits_x <- 13

#Define os limites do eixo y
SmallestP <- min(resultados$padj[which(resultados$padj > 0)], na.rm=TRUE)
axislimits_y <- ceiling(-log10(SmallestP))+0.05

# Adjust pvalues to avoid NAs
VolcanoPValues <- -log10(resultados$padj + SmallestP)
VolcanoPValues[which(is.na(VolcanoPValues))] <- 0

# Simple volcano plot
png(file=paste("./05_Results/CellLines/","VolcanoPlot_Basic_Cell_TvsTN.png", sep=""), width=2000, height=1700, res = 300)
par(mar=c(4,4,3,2), mgp=c(2,.7,0), tck=-.01)
plot(resultados$log2FoldChange, VolcanoPValues,
     xlim=c(-axislimits_x, axislimits_x), ylim=c(0, axislimits_y),
     xlab="log2 Fold Change", ylab="-log10 p-value",
     main="Volcano Plot", cex=3.5, cex.lab=2, cex.axis=2, cex.main=1, cex.sub=2,
     pch=16, col=resultados$Colour)

abline(v = log2(2), lwd=6, col=rgb(222,22,22,150,maxColorValue=255))
abline(v = -log2(2), lwd=6, col=rgb(56,50,237,150,maxColorValue=255))

text(-10, 40, labels=
       length(which(resultados$log2FoldChange < 0 & resultados$padj < 0.05)),
     cex = 1)
text(-10, 50, labels=
       length(which(resultados$log2FoldChange < -1 & resultados$padj < 0.05)),
     cex = 1)

text(10, 45, labels=
       length(which(resultados$log2FoldChange > 0 & resultados$padj < 0.05)),
     cex = 1)
text(10, 55, labels=
       length(which(resultados$log2FoldChange > 1 & resultados$padj < 0.05)),
     cex = 1)

dev.off()

# Fix Colnames
resultados$Colour <- NULL
colnames(resultados) <- paste(colnames(resultados),"_Cell",sep="")

# Create the both significant column
resultados$CellTarget <- "No"

resultados$CellTarget[
  which(resultados$padj < 0.05 & abs(resultados$log2FoldChange) > 1)
  ] <- "Yes"

#################### Naming correction TCGA data
#################### START

# Create the recipient table
Cell_naming <- matrix(ncol=6, nrow=length(rownames(resultados)))
rownames(Cell_naming) <- rownames(resultados)
colnames(Cell_naming) <- c("UCSC_symbol_split","geneID_fromSymbol", "geneID_fromAlias", "geneSymbol_fromAliasID","geneSymbol_final","geneID_final")
Cell_naming <- data.frame(Cell_naming)

#Fill with original values
Cell_naming$UCSC_symbol_split <- rownames(resultados)

# Get the GeneIdList
GeneIDList <- data.frame(toTable(org.Hs.egSYMBOL))
GeneAliasList <- data.frame(toTable(org.Hs.egALIAS2EG))

#
k <- 0
for(CellNow in rownames(Cell_naming)){
  k <- k+1
  # Use ID to find Symbol
  Position <- match(Cell_naming$UCSC_symbol_split[k], GeneIDList$symbol)
  GeneIDList$gene_id[Position] -> Cell_naming$geneID_fromSymbol[k]
  
  # Use Alias to find ID
  Position <- match(Cell_naming$UCSC_symbol_split[k], GeneAliasList$alias_symbol)
  GeneAliasList$gene_id[Position] -> Cell_naming$geneID_fromAlias[k] 
  
  # Use ID from Alias to find Symbol
  Position <- match(Cell_naming$geneID_fromAlias[k], GeneIDList$gene_id)
  GeneIDList$symbol[Position] -> Cell_naming$geneSymbol_fromAliasID[k]

  #
  if(is.na(Cell_naming$geneID_fromSymbol[k])){
    if(is.na(Cell_naming$geneID_fromAlias[k])){
      Cell_naming$UCSC_symbol_split[k] ->  Cell_naming$geneSymbol_final[k]
      Cell_naming$UCSC_symbol_split[k] ->  Cell_naming$geneID_final[k]
    } else {
      Cell_naming$geneID_fromAlias[k] ->  Cell_naming$geneID_final[k]
      Cell_naming$geneSymbol_fromAliasID[k] ->  Cell_naming$geneSymbol_final[k]
      
    }
  } else {
    Cell_naming$geneID_fromSymbol[k] ->  Cell_naming$geneID_final[k]
    Cell_naming$UCSC_symbol_split[k] ->  Cell_naming$geneSymbol_final[k]
  }
}

# Puts old names on duplicates
DuplicatedCell <- Cell_naming$geneSymbol_final[which(duplicated(Cell_naming$geneSymbol_final))]
#Solve the duplicates by the 'hard' way
for(carrot in 1:length(DuplicatedCell)){
  match(Cell_naming$geneSymbol_final, DuplicatedCell[carrot]) -> CellNowMatch
  which(!is.na(CellNowMatch)) -> NewPositions
  for(strawberry in 1:length(NewPositions)){
    Cell_naming$UCSC_symbol_split[NewPositions[strawberry]] -> Cell_naming$geneSymbol_final[NewPositions[strawberry]]
    Cell_naming$UCSC_symbol_split[NewPositions[strawberry]] -> Cell_naming$geneID_final[NewPositions[strawberry]]
  }
}

# Perform a colbind with results table
ResultsNamed <- cbind(Cell_naming[,c("geneSymbol_final", "geneID_final")], resultados)

# Write and read it again
write.table(ResultsNamed, file="./05_Results/CellLines/results_CellLines_TvsNT_results_DESeq2.txt")
write.csv2(ResultsNamed, file="./05_Results/CellLines/results_CellLines_TvsNT_results_DESeq2.csv")
ResultsNamed <- data.frame(read.table(file="./05_Results/CellLines/results_CellLines_TvsNT_results_DESeq2.txt", stringsAsFactors = FALSE))


######### PATHWAY ENRICHMENT

#################### GO ENRICHMENT
#################### START

# Get DE genes in both comparisons
ResultsNamed_DE <- ResultsNamed[which(ResultsNamed$CellTarget == "Yes"), ]

# Create the Gene Vector
# Gene vector is a simple vector stating 1 for DE genes and
# 0 for non-de genes.
gene.vector_DE <- as.integer(ResultsNamed$geneSymbol_final%in%ResultsNamed_DE$geneSymbol_final)

# As this vector is unamed, put the right names on it
names(gene.vector_DE) <- ResultsNamed$geneSymbol_final

#Probability Weighting Function for Gene length
# Also plot the model
# the second position is for genome version and the third is for ID type
# hg19 is TCGA assembly version
# geneSymbol is intended to Gene Names
# knownGene is intended to GeneIDs
# To check all human genomes and supported ids, uncomment:
# supportedGenomes()[which(supportedGenomes()[,"species"] == "Human"),]
png(file=paste("./05_Results/CellLines/Pathway/GO/ProbabilityWeightingFunction",".png", sep=""), width=2000, height=2000, res = 400)
pwf <- nullp(gene.vector_DE,"hg19","geneSymbol")
dev.off()

# Perform Wallenius GO Analysis
GO.wall <- goseq(pwf, "hg19", "geneSymbol", 
                 test.cats=c("GO:CC", "GO:BP", "GO:MF"),
                 method = "Wallenius", use_genes_without_cat=FALSE)

# Calculate FDR for Over represented GOs (the most important value)
GO.wall$over_represented_fdr <- p.adjust(GO.wall$over_represented_pvalue, method="BH")

# Calculate FDR for Under represented GOs
GO.wall$under_represented_fdr <- p.adjust(GO.wall$under_represented_pvalue, method="BH")

# Split GO table for each Go ontology
GO.wall_CC <- GO.wall[which(GO.wall[,"ontology"]=='CC'),]
GO.wall_BP <- GO.wall[which(GO.wall[,"ontology"]=='BP'),]
GO.wall_MF <- GO.wall[which(GO.wall[,"ontology"]=='MF'),]

# Write the files
write.table(GO.wall_CC, file=paste("./05_Results/CellLines/Pathway/GO/CC",".txt",sep=""), sep="\t")
write.table(GO.wall_BP, file=paste("./05_Results/CellLines/Pathway/GO/BP",".txt",sep=""), sep="\t")
write.table(GO.wall_MF, file=paste("./05_Results/CellLines/Pathway/GO/MF",".txt",sep=""), sep="\t")

# Split GO table for each Go ontology and FDR values
GO.wall_fdr_CC <- GO.wall[which(GO.wall[,"ontology"]=='CC' & GO.wall[,"over_represented_fdr"] < 0.05),]
GO.wall_fdr_BP <- GO.wall[which(GO.wall[,"ontology"]=='BP' & GO.wall[,"over_represented_fdr"] < 0.05),]
GO.wall_fdr_MF <- GO.wall[which(GO.wall[,"ontology"]=='MF' & GO.wall[,"over_represented_fdr"] < 0.05),]

# Define category order
MainNames <- c("Cellular component", "Biological Process", "Molecular Function")

# Plot First 15 cat in each
i <- 0
for(pineapple in c("GO.wall_fdr_CC", "GO.wall_fdr_BP", "GO.wall_fdr_MF")){
  i <- i+1
  get(pineapple) -> PlotNowGO
  
  png(file=paste("./05_Results/CellLines/Pathway/GO/GOEnrichPlot_15first_",pineapple,".png",sep=""),
      width=6000, height=(660+(110*length(PlotNowGO$ontology[1:15]))), res = 400)
  par(mar=c(5,45,2,2),lwd = 2.5)
  
  barplot(-log10(as.numeric(PlotNowGO[1:15,"over_represented_fdr"])), 
          names=PlotNowGO[1:15,"term"], las=2, horiz=TRUE, xlab="-log(FDR)",
          main=MainNames[i], lwd=2, cex.lab=1.7, cex.axis=2, cex.main=2,
          axes=FALSE, col=brewer.pal(8,"Set1")[i], 
          border="white",
          cex.names=1.5, space=0.001)
  
  abline(v=1:floor(-log10(min(PlotNowGO[,"over_represented_fdr"]))),
         col="white", lwd=4)
  
  axis(side = 1, lwd = 1.5, cex.axis=1.4)
  
  dev.off()
}

# Plot All data
i <- 0
for(pineapple in c("GO.wall_fdr_CC", "GO.wall_fdr_BP", "GO.wall_fdr_MF")){
  i <- i+1
  get(pineapple) -> PlotNowGO
  
  png(file=paste("./05_Results/CellLines/Pathway/GO/GOEnrichPlot_all_",pineapple,".png",sep=""),
      width=5000, height=(660+(110*length(PlotNowGO$ontology))), res = 400)
  par(mar=c(5,40,2,2),lwd = 2.5)
  
  barplot(-log10(as.numeric(PlotNowGO[,"over_represented_fdr"])), 
          names=PlotNowGO[,"term"], las=2, horiz=TRUE, xlab="-log(FDR)",
          main=MainNames[i], lwd=2, cex.lab=1.7, cex.axis=2, cex.main=2,
          axes=FALSE, col=brewer.pal(8,"Set1")[i], 
          border="white",
          cex.names=1.5, space=0.001)
  
  abline(v=1:floor(-log10(min(PlotNowGO[,"over_represented_fdr"]))),
         col="white", lwd=4)
  
  axis(side = 1, lwd = 1.5, cex.axis=1.4)
  
  dev.off()
}

#################### GO ENRICHMENT
#################### END

#################### REACTOME ENRICHMENT
#################### START


# Another tool for conversion from gene name to gene ID, after HUGO translation
geneIDList_DE <- bitr(ResultsNamed_DE$geneSymbol_final, fromType="SYMBOL",
                      toType="ENTREZID", OrgDb="org.Hs.eg.db")
geneIDList_Universe <- bitr(ResultsNamed$geneSymbol_final, fromType="SYMBOL",
                            toType="ENTREZID", OrgDb="org.Hs.eg.db")

# perform enrichment in REACTOME
REACTOME_Enriched <- enrichPathway(geneIDList_DE$ENTREZID,
                                   organism = "human",
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 1,
                                   universe = geneIDList_Universe$ENTREZID,
                                   minGSSize = 2,
                                   readable = TRUE)

# Makes summary
REACTOME_Enriched_summary <- summary(REACTOME_Enriched)

# Write
write.table(REACTOME_Enriched_summary, file=paste("./05_Results/CellLines/Pathway/REACTOME/REAC.enrichment",".txt",sep=""), sep="\t")

#Get fdr <0.05
REACTOME_Enriched_fdr <- REACTOME_Enriched_summary[which(REACTOME_Enriched_summary$p.adjust<0.05),]


# REACTOME Plot
# First 15
REACTOME_Enriched_fdr[1:15,]-> PlotNowREACT

png(file=paste("./05_Results/CellLines/Pathway/REACTOME/REACT_EnrichPlot_15first.png",sep=""),
    width=6000, height=(660+(110*length(PlotNowREACT$Description))), res = 450)
par(mar=c(5,45,2,2),lwd = 2.5)

barplot(-log10(as.numeric(PlotNowREACT[,"p.adjust"])), 
        names=PlotNowREACT[,"Description"], las=2, horiz=TRUE, xlab="-log(FDR)",
        main="REACTOME", lwd=2, cex.lab=1.7, cex.axis=2, cex.main=2,
        axes=FALSE, col=brewer.pal(8,"Set1")[4], 
        border="white",
        cex.names=1.5, space=0.001)

abline(v=1:floor(-log10(min(PlotNowREACT[,"p.adjust"]))),
       col="white", lwd=4)

axis(side = 1, lwd = 1.5, cex.axis=1.4)

dev.off()

#All
REACTOME_Enriched_fdr[,]-> PlotNowREACT

png(file=paste("./05_Results/CellLines/Pathway/REACTOME/REACT_EnrichPlot_All.png",sep=""),
    width=6000, height=(660+(110*length(PlotNowREACT$Description))), res = 450)
par(mar=c(5,45,2,2),lwd = 2.5)

barplot(-log10(as.numeric(PlotNowREACT[,"p.adjust"])), 
        names=PlotNowREACT[,"Description"], las=2, horiz=TRUE, xlab="-log(FDR)",
        main="REACTOME", lwd=2, cex.lab=1.7, cex.axis=2, cex.main=2,
        axes=FALSE, col=brewer.pal(8,"Set1")[4], 
        border="white",
        cex.names=1.5, space=0.001)

abline(v=1:floor(-log10(min(PlotNowREACT[,"p.adjust"]))),
       col="white", lwd=4)

axis(side = 1, lwd = 1.5, cex.axis=1.4)

dev.off()

# perform enrichment in REACTOME with cutoffs
REACTOME_Enriched_2 <- enrichPathway(geneIDList_DE$ENTREZID,
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 1,
                                     universe = geneIDList_Universe$ENTREZID,
                                     minGSSize = 2,
                                     readable = TRUE)

# Plot Enrichment Map without names
png(file=paste("./05_Results/CellLines/Pathway/REACTOME/REACT_EnrichMap_NoNames.png",sep=""),
    width=6000, heigh=6000, res=1)
enrichMap(REACTOME_Enriched,
          layout=igraph::layout.kamada.kawai,
          vertex.label.cex = 1,
          vertex.label.font = 2, 
          #vertex.label.color = "#666666",
          foldChange = NULL,
          fixed = TRUE,
          col.bin = 10)
dev.off()

# Plot Enrichment Map With Names
png(file=paste("./05_Results/CellLines/Pathway/REACTOME/REACT_EnrichMap.png",sep=""),
    width=6000, heigh=6000, res=450)
enrichMap(REACTOME_Enriched,
          layout=igraph::layout.kamada.kawai,
          vertex.label.cex = 1,
          vertex.label.font = 2, 
          #vertex.label.color = "#666666",
          foldChange = NULL,
          fixed = TRUE,
          col.bin = 10)
dev.off()

# Plot the cnetplot for the First Cat
png(file=paste("./05_Results/CellLines/Pathway/REACTOME/REACT_cnetplot.png",sep=""),
    width=6000, heigh=6000, res=450)
cnetplot(REACTOME_Enriched, categorySize="Count", foldChange=NULL,
         showCategory = 1)
dev.off()


#################### REACTOME
#################### END

#################### Disease Ontology Semantic and Enrichment analysis
#################### START

# Another tool for conversion from gene name to gene ID, after HUGO translation
geneIDList_DE <- bitr(ResultsNamed_DE$geneSymbol_final, fromType="SYMBOL",
                      toType="ENTREZID", OrgDb="org.Hs.eg.db")
geneIDList_Universe <- bitr(ResultsNamed$geneSymbol_final, fromType="SYMBOL",
                            toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Perform DO Enrichment analysis
DO_Enrichment <- enrichDO(gene = geneIDList_DE$ENTREZID,
                          ont           = "DO", 
                          pvalueCutoff  = 1,
                          pAdjustMethod = "BH",
                          universe      = geneIDList_Universe$ENTREZID, 
                          minGSSize     = 2,
                          qvalueCutoff  = 1,
                          readable      = TRUE)

# Makes summary
DO_Enriched_summary <- summary(DO_Enrichment)

# Write
write.table(DO_Enriched_summary, file=paste("./05_Results/CellLines/Pathway/DOSE/DO.enrichment",".txt",sep=""), sep="\t")

#Get fdr <0.05
DO_Enriched_fdr <- DO_Enriched_summary[which(DO_Enriched_summary$p.adjust<0.05),]


# DO Plot
# First 15
DO_Enriched_fdr[1:15,]-> PlotNowDO

png(file=paste("./05_Results/CellLines/Pathway/DOSE/DO_EnrichPlot_15first.png",sep=""),
    width=6000, height=(660+(110*length(PlotNowDO$Description))), res = 450)
par(mar=c(5,45,2,2),lwd = 2.5)

barplot(-log10(as.numeric(PlotNowDO[,"p.adjust"])), 
        names=PlotNowDO[,"Description"], las=2, horiz=TRUE, xlab="-log(FDR)",
        main="Disease Ontology", lwd=2, cex.lab=1.7, cex.axis=2, cex.main=2,
        axes=FALSE, col=brewer.pal(8,"Set1")[5], 
        border="white",
        cex.names=1.5, space=0.001)

abline(v=1:floor(-log10(min(PlotNowDO[,"p.adjust"]))),
       col="white", lwd=4)
axis(side = 1, lwd = 1.5, cex.axis=1.4)
dev.off()

# DO Plot
# All
DO_Enriched_fdr[,]-> PlotNowDO

png(file=paste("./05_Results/CellLines/Pathway/DOSE/DO_EnrichPlot_All.png",sep=""),
    width=6000, height=(660+(110*length(PlotNowDO$Description))), res = 450)
par(mar=c(5,45,2,2),lwd = 2.5)

barplot(-log10(as.numeric(PlotNowDO[,"p.adjust"])), 
        names=PlotNowDO[,"Description"], las=2, horiz=TRUE, xlab="-log(FDR)",
        main="Disease Ontology", lwd=2, cex.lab=1.7, cex.axis=2, cex.main=2,
        axes=FALSE, col=brewer.pal(8,"Set1")[5], 
        border="white",
        cex.names=1.5, space=0.001)

abline(v=1:floor(-log10(min(PlotNowDO[,"p.adjust"]))),
       col="white", lwd=4)

axis(side = 1, lwd = 1.5, cex.axis=1.4)

dev.off()

#################### Disease Ontology Semantic and Enrichment analysis
#################### END

#################### PathWay Enrichment Cross
#################### START

i <- 0
for(TypeNow in c("CC","BP","MF")){

i <- i+1
# Read the Tables
FromTCGA <- data.frame(read.table(file=paste("./04_Results/GO/",TypeNow,".txt",sep=""), stringsAsFactors = FALSE))
FromCELL <- data.frame(read.table(file=paste("./05_Results/CellLines/Pathway/GO/",TypeNow,".txt",sep=""), stringsAsFactors = FALSE))

# Intersect the significant values
FromTCGA_FDR<- FromTCGA$category[which(FromTCGA$over_represented_fdr < 0.05)]
FromCELL_FDR <- FromCELL$category[which(FromCELL$over_represented_fdr < 0.05)]

# Select significant and intersected
FromTCGA <- FromTCGA[match(intersect(FromTCGA_FDR, FromCELL_FDR), FromTCGA$category),]
FromCELL <- FromCELL[match(intersect(FromTCGA_FDR, FromCELL_FDR), FromCELL$category),]

png(file=paste("./05_Results/Intersection/Pathway/GOEnrichPlot_all_",TypeNow,"_TCGA.png",sep=""),
    width=5000, height=(660+(110*length(FromTCGA$ontology))), res = 400)
par(mar=c(5,40,2,2),lwd = 2.5)

barplot(-log10(as.numeric(FromTCGA[,"over_represented_fdr"])), 
        names=FromTCGA[,"term"], las=2, horiz=TRUE, xlab="-log(FDR)",
        main=MainNames[i], lwd=2, cex.lab=1.7, cex.axis=2, cex.main=2,
        axes=FALSE, col=brewer.pal(8,"Set1")[i], 
        border="white",
        cex.names=1.5, space=0.001)

abline(v=1:floor(-log10(min(FromTCGA[,"over_represented_fdr"]))),
       col="white", lwd=4)

axis(side = 1, lwd = 1.5, cex.axis=1.4)

dev.off()

png(file=paste("./05_Results/Intersection/Pathway/GOEnrichPlot_all_",TypeNow,"_Cell.png",sep=""),
    width=5000, height=(660+(110*length(FromTCGA$ontology))), res = 400)
par(mar=c(5,40,2,2),lwd = 2.5)

barplot(-log10(as.numeric(FromCELL[,"over_represented_fdr"])), 
        names=FromCELL[,"term"], las=2, horiz=TRUE, xlab="-log(FDR)",
        main=MainNames[i], lwd=2, cex.lab=1.7, cex.axis=2, cex.main=2,
        axes=FALSE, col=brewer.pal(8,"Set1")[i], 
        border="white",
        cex.names=1.5, space=0.001)

abline(v=1:floor(-log10(min(FromCELL[,"over_represented_fdr"]))),
       col="white", lwd=4)

axis(side = 1, lwd = 1.5, cex.axis=1.4)
dev.off()
}

#################### PathWay Enrichment Cross
#################### END


# Boxplot genes
# This should take long
# Set the levels
mylevels <- c("NonTriple", "Triple")

##################################### >>>>> ONLY DE GENES
# Create table just with DE and with FC cutoff
resultadosDE <- ResultsNamed[which(ResultsNamed$CellTarget == "Yes"),]

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
  PlottingNow <- data.frame(matrix(nrow=dim(NormalizedGeneExpression)[2], ncol=0))
  rownames(PlottingNow) <- colnames(NormalizedGeneExpression)
  PlottingNow$value <- as.numeric(log2(NormalizedGeneExpression[Target,]+1))
  PlottingNow$condition <- colData$conditions
  
  # perform the boxplot itself
  png(paste("./05_Results/CellLines/Boxplots_DEGenes/Cell_",alvos,".png", sep=""), width=5000, height=7000, res=1200, type='cairo')
  boxplot(PlottingNow$value~PlottingNow$condition,
          las = 1,
          col = c("#80D99A","#79D8DB"),
          border = c("#3AC664","#00BFC4"),
          outline=TRUE,
          range=0.5,
          boxwex=0.8,
          notch = FALSE,
          lwd=2,
          #width=levelProportions,
          ylab ="log2(RSEM+1)",
          main=alvos,
          par(cex.lab=1.5, cex.axis=1),
          outpch=NA
  )

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

# Boxplot ALL genes
# This should take long
# Set the levels
mylevels <- c("NonTriple", "Triple")

# Create table just with DE and with FC cutoff
resultadosDE <- ResultsNamed

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
  PlottingNow <- data.frame(matrix(nrow=dim(NormalizedGeneExpression)[2], ncol=0))
  rownames(PlottingNow) <- colnames(NormalizedGeneExpression)
  PlottingNow$value <- as.numeric(log2(NormalizedGeneExpression[Target,]+1))
  PlottingNow$condition <- colData$conditions
  
  # perform the boxplot itself
  png(paste("./05_Results/CellLines/Boxplots_AllGenes/Cell_",alvos,".png", sep=""), width=4000, height=7000, res=1200, type='cairo')
  boxplot(PlottingNow$value~PlottingNow$condition,
          las = 1,
          col = c("#80D99A","#79D8DB"),
          border = c("#3AC664","#00BFC4"),
          outline=TRUE,
          range=0.5,
          boxwex=0.8,
          notch = FALSE,
          lwd=2,
          #width=levelProportions,
          ylab ="log2(RSEM+1)",
          main=alvos,
          par(cex.lab=1.5, cex.axis=1),
          outpch=NA
  )

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


########################
# qPCR Validation

# some definitions bellow
qPCR <- read.table(file="./ExternalFiles/qPCR_Validation.txt", header=TRUE, stringsAsFactors = FALSE, sep ="\t")
RNASeqOrder <- c("BT549_LNBio","MCF7_LNBio","MDAMB231_LNBio","MDAMB436_LNBio","MDAMB468_LNBio","SKBR3_LNBio")
qPCROrder <- c("BT549","MCF7","MDAMB231","MDAMB436","MDAMB468","SKBR3")
OutlierRemovalVector <- paste(RNASeqOrder, "_Genes",sep="")

# Pearson Correlations
for(cabbage in 1:length(RNASeqOrder)){
  
  #
  qPCR[which(qPCR$Cell == qPCROrder[cabbage]),] -> qPCR_Now
  
  # Fills with zeros
  qPCR_Now$RNAseq <- 0
  
  # Loop for table filling
  for(strawberry in 1:length(rownames(qPCR_Now))){
    qPCR_Now$RNAseq[strawberry] <- NormalizedGeneExpression[qPCR_Now$Target[strawberry],RNASeqOrder[cabbage]] 
  }
  
  # Plotting
  1/qPCR_Now$dCt -> x
  log2(qPCR_Now$RNAseq+1) -> y
  w <- round(cor(x,y),4)

  png(file=paste("./05_Results/CellLines/qPCR/LogRSEM_Inverse-dCt_", RNASeqOrder[cabbage] ,".png", sep=""), res=350, height = 1500, width = 2000)
  plot(x,y,
       ylab = "log2(RSEM+1) - RNA-seq",
       xlab = "1/??Ct - qPCR",
       cex=1.5, cex.axis=1.4,
       cex.lab=1.5, cex.main=1.6,
       col=rgb(50,50,50,50,maxColorValue=255),
       xlim=c(0,0.1),
       ylim=c(0,16),
       main=RNASeqOrder[cabbage],
       pch=16)
  abline(lm(y~x), 
         col="red", lwd=2, lty=4)
  text(0.08, 1, paste("Pearson = ",w,sep=""), cex=1.2, font=3)
  text(0.08, 3.5, paste("r-squared = ",
                    round(summary(lm(y~x))$adj.r.squared,4),
                    sep=""), cex=1.2, font=3)
  dev.off()
  
  
  # Plotting without small RSEMs
  qPCR_Now[qPCR_Now$RNAseq>50,] -> qPCR_Now
  
  1/qPCR_Now$dCt -> x
  log2(qPCR_Now$RNAseq+1) -> y
  w <- round(cor(x,y),4)

  png(file=paste("./05_Results/CellLines/qPCR/LogRSEM_Inverse-dCt_RSEM_50_", RNASeqOrder[cabbage] ,".png", sep=""), res=350, height = 1500, width = 2000)
  plot(x,y,
       ylab = "log2(RSEM+1) - RNA-seq",
       xlab = "1/??Ct - qPCR",
       cex=1.5, cex.axis=1.4,
       cex.lab=1.5, cex.main=1.6,
       col=rgb(50,50,50,50,maxColorValue=255),
       xlim=c(0,0.1),
       ylim=c(0,16),
       main=RNASeqOrder[cabbage],
       pch=16)
  abline(lm(y~x), 
         col="red", lwd=2, lty=4)
  text(0.08, 1, paste("Pearson = ",w,sep=""), cex=1.2, font=3)
  text(0.08, 3.5, paste("r-squared = ",
                        round(summary(lm(y~x))$adj.r.squared,4),
                        sep=""), cex=1.2, font=3)
  dev.off()
  
  
  # Outlier removal
  assign(OutlierRemovalVector[cabbage],
         unique(qPCR_Now[which(y/x > mean(y/x)-sd(y/x) & y/x < mean(y/x)+sd(y/x)),]$Target))
}


# Spearman Correlations
for(cabbage in 1:length(RNASeqOrder)){
  
  #
  qPCR[which(qPCR$Cell == qPCROrder[cabbage]),] -> qPCR_Now
  
  # Fills with zeros
  qPCR_Now$RNAseq <- 0
  
  # Loop for table filling
  for(strawberry in 1:length(rownames(qPCR_Now))){
    qPCR_Now$RNAseq[strawberry] <- NormalizedGeneExpression[qPCR_Now$Target[strawberry],RNASeqOrder[cabbage]] 
  }
  
  # Plotting
  1/qPCR_Now$dCt -> x
  log2(qPCR_Now$RNAseq+1) -> y
  w <- round(cor(x,y, method="spearman"),4)
  
  png(file=paste("./05_Results/CellLines/qPCR_Spearman/LogRSEM_Inverse-dCt_", RNASeqOrder[cabbage] ,".png", sep=""), res=350, height = 1500, width = 2000)
  plot(x,y,
       ylab = "log2(RSEM+1) - RNA-seq",
       xlab = "1/??Ct - qPCR",
       cex=1.5, cex.axis=1.4,
       cex.lab=1.5, cex.main=1.6,
       col=rgb(50,50,50,50,maxColorValue=255),
       xlim=c(0,0.1),
       ylim=c(0,16),
       main=RNASeqOrder[cabbage],
       pch=16)
  abline(lm(y~x), 
         col="red", lwd=2, lty=4)
  text(0.08, 1, paste("Spearman = ",w,sep=""), cex=1.2, font=3)
  text(0.08, 3.5, paste("r-squared = ",
                        round(summary(lm(y~x))$adj.r.squared,4),
                        sep=""), cex=1.2, font=3)
  dev.off()
  
}   


########## RSEM DENSITIES PLOT
cols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(colnames(GeneExpression)))
data.frame(GeneExpression) -> GeneExpression
data.frame(NormalizedGeneExpression) -> NormalizedGeneExpression

#Plot Normalized RSEM densities
png(file="./05_Results/CellLines/Densities/RSEM_Cells_Norm.png", width=4000, height = 3000, res=600)
# Start the plot itself by the first column
plot(density(log2(NormalizedGeneExpression[,1])),
     main= "Normalized log2(RSEM) Density",
     ylab="Density",
     xlab="",
     col=cols[1],
     lwd=3,
     cex=1.5,
     cex.lab=1.5,
     yaxt="n",
     xaxt="n",
     bty="n",
     ylim=c(0,0.09)
)

#Add all other columns, one by one
for(ActualPatient in 2:length(colnames(GeneExpression))){
  lines(density(log2(NormalizedGeneExpression[,ActualPatient])),
        col=cols[ActualPatient],
        lwd=3)  
}

# Add the axis
axis(side = 1, lwd = 2,cex.axis=1.4)
axis(side = 2, lwd = 2,cex.axis=1.4)
box(lwd=2)

dev.off()

#Plot Raw RSEM densities
png(file="./05_Results/CellLines/Densities/RSEM_Cells_Raw.png", width=4000, height = 3000, res=600)
# Start the plot itself by the first column
plot(density(log2(GeneExpression[,1])),
     main= "log2(RSEM) Density",
     ylab="Density",
     xlab="",
     col=cols[1],
     lwd=3,
     cex=1.5,
     cex.lab=1.5,
     yaxt="n",
     xaxt="n",
     bty="n",
     ylim=c(0,0.09)
)

#Add all other columns, one by one
for(ActualPatient in 2:length(colnames(GeneExpression))){
  lines(density(log2(GeneExpression[,ActualPatient])),
        col=cols[ActualPatient],
        lwd=3)  
}

# Add the axis
axis(side = 1, lwd = 2,cex.axis=1.4)
axis(side = 2, lwd = 2,cex.axis=1.4)
box(lwd=2)

dev.off()


################## RSEM PLOTS JUST SOME CELL LINES
cols <- brewer.pal(8,"Set1")
png(file="./05_Results/CellLines/Densities/RSEM_LNBioCells_Raw.png", width=4000, height = 3000, res=600)
plot(density(log2(GeneExpression$BT549_LNBio)),
     main= "Raw log2(RSEM+1) Density",
     ylab="Density",
     xlab="",
     col=cols[1],
     lwd=3,
     cex=1.5,
     cex.lab=1.5,
     yaxt="n",
     xaxt="n",
     ylim=c(0,0.09)
)

lines(density(log2(GeneExpression$MCF7_LNBio)),
      col=cols[2],
      lwd=3)

lines(density(log2(GeneExpression$MDAMB231_LNBio)),
      col=cols[3],
      lwd=3)


lines(density(log2(GeneExpression$MDAMB468_LNBio)),
      col=cols[4],
      lwd=3)

lines(density(log2(GeneExpression$MDAMB436_LNBio)),
      col=cols[5],
      lwd=3)

lines(density(log2(GeneExpression$SKBR3_LNBio)),
      col=cols[6],
      lwd=3)


legend(-7,0.094,
       c("BT549", "MCF7","MDAMB436","MDAMB468","MDAMB231","SKBR3"),
       lty=c(1,1,1,1,1,1),
       lwd=c(3,3,3,3,3,3),
       col=cols[1:6])

axis(side = 1, lwd = 2,cex.axis=1.4)
axis(side = 2, lwd = 2,cex.axis=1.4)
box(lwd=2)

dev.off()


png(file="./05_Results/CellLines/Densities/RSEM_LNBioCells_Norm.png", width=4000, height = 3000, res=600)
plot(density(log2(NormalizedGeneExpression$BT549_LNBio)),
     main= "Normalized log2(RSEM) Density",
     ylab="Density",
     xlab="",
     col=cols[1],
     lwd=3,
     cex=1.5,
     cex.lab=1.5,
     yaxt="n",
     xaxt="n",
     ylim=c(0,0.095)
)

lines(density(log2(NormalizedGeneExpression$MCF7_LNBio)),
      col=cols[2],
      lwd=3)

lines(density(log2(NormalizedGeneExpression$MDAMB231_LNBio)),
      col=cols[3],
      lwd=3)

lines(density(log2(NormalizedGeneExpression$MDAMB468_LNBio)),
      col=cols[4],
      lwd=3)

lines(density(log2(NormalizedGeneExpression$MDAMB436_LNBio)),
      col=cols[5],
      lwd=3)

lines(density(log2(NormalizedGeneExpression$SKBR3_LNBio)),
      col=cols[6],
      lwd=3)

legend(-7,0.094,
       c("BT549", "MCF7","MDAMB436","MDAMB468","MDAMB231","SKBR3"),
       lty=c(1,1,1,1,1,1),
       lwd=c(3,3,3,3,3,3),
       col=cols[1:6])

axis(side = 1, lwd = 2,cex.axis=1.4)
axis(side = 2, lwd = 2,cex.axis=1.4)
box(lwd=2)

dev.off()

############ CROSSING LISTS START

# Read back expression table from TCGA
ResultsCell <- data.frame(read.table(file="./05_Results/CellLines/results_CellLines_TvsNT_results_DESeq2.txt", stringsAsFactors = FALSE))
ResultsTCGA <- data.frame(read.table(file="./03_Results/results_DESeq2_intersected_naming.txt", stringsAsFactors = FALSE))

# Get the intersection gene list
IntersectionGenes <- intersect(ResultsCell$geneSymbol_final,
ResultsTCGA$geneSymbol_final)

# Joint tables in right order
ResultsMerged <- cbind(
ResultsTCGA[match(IntersectionGenes, ResultsTCGA$geneSymbol_final),],
ResultsCell[match(IntersectionGenes, ResultsCell$geneSymbol_final),])

# Remove redundant and old columns
ResultsMerged[,17] <- NULL
ResultsMerged[,16] <- NULL
ResultsMerged$BothTarget -> ResultsMerged$TCGATarget
ResultsMerged$BothTarget <- NULL

# Create column for final target assignment
ResultsMerged$TCGA_Cell_Target <- "Irrelevant"

# Up regulated in everyone and significant
ResultsMerged$TCGA_Cell_Target[
  which(ResultsMerged$CellTarget == "Yes" &
          ResultsMerged$TCGATarget == "Both" &
          ResultsMerged$log2FoldChange_TvsNT > 1 &
          ResultsMerged$log2FoldChange_TvsH > 1 &
          ResultsMerged$log2FoldChange_Cell > 1)
  ] <- "All_Upregulated"

# Down regulated in everyone and significant
ResultsMerged$TCGA_Cell_Target[
  which(ResultsMerged$CellTarget == "Yes" &
          ResultsMerged$TCGATarget == "Both" &
          ResultsMerged$log2FoldChange_TvsNT < -1 &
          ResultsMerged$log2FoldChange_TvsH < -1 &
          ResultsMerged$log2FoldChange_Cell < -1)
  ] <- "All_Downregulated"

#
ResultsMerged$TCGA_Cell_Target[
which(!ResultsMerged$TCGATarget == "Both" &
        !ResultsMerged$TCGATarget == "vsHealthy" &
        !ResultsMerged$TCGATarget == "vsNonTriple" &
        ResultsMerged$CellTarget == "Yes" &
        ResultsMerged$log2FoldChange_Cell > 1)
] <- "CellOnly_Upregulated"

#
ResultsMerged$TCGA_Cell_Target[
  which(!ResultsMerged$TCGATarget == "Both" &
          !ResultsMerged$TCGATarget == "vsHealthy" &
          !ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "Yes" &
          ResultsMerged$log2FoldChange_Cell < -1)
  ] <- "CellOnly_Downregulated"

ResultsMerged$TCGA_Cell_Target[
  which(!ResultsMerged$TCGATarget == "Both" &
          ResultsMerged$TCGATarget == "vsHealthy" &
          !ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "Yes" &
          ResultsMerged$log2FoldChange_Cell > 1 &
          ResultsMerged$log2FoldChange_TvsH > 1)
  ] <- "Cell_Healthy_Upregulated"

ResultsMerged$TCGA_Cell_Target[
  which(!ResultsMerged$TCGATarget == "Both" &
          ResultsMerged$TCGATarget == "vsHealthy" &
          !ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "Yes" &
          ResultsMerged$log2FoldChange_Cell < -1 &
          ResultsMerged$log2FoldChange_TvsH < -1)
  ] <- "Cell_Healthy_Downregulated"

ResultsMerged$TCGA_Cell_Target[
  which(!ResultsMerged$TCGATarget == "Both" &
          !ResultsMerged$TCGATarget == "vsHealthy" &
          ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "Yes" &
          ResultsMerged$log2FoldChange_Cell > 1 &
        ResultsMerged$log2FoldChange_TvsNT > 1)
  ] <- "Cell_TCGA_TN_Upregulated"

ResultsMerged$TCGA_Cell_Target[
  which(!ResultsMerged$TCGATarget == "Both" &
          !ResultsMerged$TCGATarget == "vsHealthy" &
          ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "Yes" &
          ResultsMerged$log2FoldChange_Cell < -1 &
        ResultsMerged$log2FoldChange_TvsNT < -1)
  ] <- "Cell_TCGA_TN_Downregulated"

ResultsMerged$TCGA_Cell_Target[
  which(ResultsMerged$TCGATarget == "Both" &
          !ResultsMerged$TCGATarget == "vsHealthy" &
          !ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "No" &
          ResultsMerged$log2FoldChange_TvsNT > 1 &
          ResultsMerged$log2FoldChange_TvsH > 1)
  ] <- "TCGAOnly_Upregulated"

ResultsMerged$TCGA_Cell_Target[
  which(ResultsMerged$TCGATarget == "Both" &
          !ResultsMerged$TCGATarget == "vsHealthy" &
          !ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "No" &
          ResultsMerged$log2FoldChange_TvsNT < -1 &
          ResultsMerged$log2FoldChange_TvsH < -1)
  ] <- "TCGAOnly_Downregulated"

ResultsMerged$TCGA_Cell_Target[
  which(!ResultsMerged$TCGATarget == "Both" &
          ResultsMerged$TCGATarget == "vsHealthy" &
          !ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "No" &
          ResultsMerged$log2FoldChange_TvsH > 1)
  ] <- "HealthyOnly_Upregulated"

ResultsMerged$TCGA_Cell_Target[
  which(!ResultsMerged$TCGATarget == "Both" &
          ResultsMerged$TCGATarget == "vsHealthy" &
          !ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "No" &
          ResultsMerged$log2FoldChange_TvsH < -1)
  ] <- "HealthyOnly_Downregulated"

ResultsMerged$TCGA_Cell_Target[
  which(!ResultsMerged$TCGATarget == "Both" &
          !ResultsMerged$TCGATarget == "vsHealthy" &
          ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "No" &
          ResultsMerged$log2FoldChange_TvsNT > 1)
  ] <- "TCGATNOnly_Upregulated"

ResultsMerged$TCGA_Cell_Target[
  which(!ResultsMerged$TCGATarget == "Both" &
          ResultsMerged$TCGATarget == "vsHealthy" &
          !ResultsMerged$TCGATarget == "vsNonTriple" &
          ResultsMerged$CellTarget == "No" &
          ResultsMerged$log2FoldChange_TvsNT < -1)
  ] <- "TCGATNOnly_Downregulated"

# Save and Read
write.table(ResultsMerged, file="./05_Results/Intersection/IntersectionAllTable.txt")
write.csv2(ResultsMerged, file="./05_Results/Intersection/IntersectionAllTable.csv")
ResultsMerged <- data.frame(read.table(file="./05_Results/Intersection/IntersectionAllTable.txt", stringsAsFactors = FALSE))

# Venn diagram
# Create lists with genes for an 7 class venn
# upregulated list
UpregulatedGenes <- list(
Cell = rownames(ResultsMerged)[
  which(ResultsMerged$padj_Cell < 0.05 &
      ResultsMerged$log2FoldChange_Cell > 1)],
Healthy = rownames(ResultsMerged)[
  which(ResultsMerged$padj_TvsH < 0.05 &
          ResultsMerged$log2FoldChange_TvsH > 1)],
TvsNT = rownames(ResultsMerged)[
  which(ResultsMerged$padj_TvsNT < 0.05 &
          ResultsMerged$log2FoldChange_TvsNT > 1)])
# Downregulated list
DownregulatedGenes <- list(
  Cell = rownames(ResultsMerged)[
    which(ResultsMerged$padj_Cell < 0.05 &
            ResultsMerged$log2FoldChange_Cell < -1)],
  Healthy = rownames(ResultsMerged)[
    which(ResultsMerged$padj_TvsH < 0.05 &
            ResultsMerged$log2FoldChange_TvsH < -1)],
  TvsNT = rownames(ResultsMerged)[
    which(ResultsMerged$padj_TvsNT < 0.05 &
            ResultsMerged$log2FoldChange_TvsNT < -1)])

# UP Regulated
venn.plot <- venn.diagram(
  UpregulatedGenes,
  height = 3000, width = 3000, resolution = 500, imagetype = "png",
  main="VennDiagramUP",
  main.cex=1,
  print.mode="raw",
  lwd=c("4","4","4"),
  col=c("#F8766D"),
  fill=c("#F8766D"),
  alpha=0.5,
  cat.cex = 1.5,
  cex=c(2,2,2,2,2,2,2),
  margin = 0.12,
  cat.dist=c(0.06,0.06,0.06),
  category.names = c("Cell","TvsH","TvsNT"),
  filename="./05_Results/Intersection/VennDiagrams/VennAll_Up.png")

# Down Regulated
venn.plot <- venn.diagram(
  DownregulatedGenes,
  height = 3000, width = 3000, resolution = 500, imagetype = "png",
  main="VennDiagramDown",
  main.cex=1,
  print.mode="raw",
  lwd=c("4","4","4"),
  col=c("#00BFC4"),
  fill=c("#00BFC4"),
  alpha=0.5,
  cat.cex = 1.5,
  cex=c(2,2,2,2,2,2,2),
  margin = 0.12,
  cat.dist=c(0.06,0.06,0.06),
  category.names = c("Cell","TvsH","TvsNT"),
  filename="./05_Results/Intersection/VennDiagrams/VennAll_Down.png")

# Correlation plots
#Generate the vector with the highest padj for the correlation
Maximum_padj <- rep(1, dim(ResultsMerged)[1])
for(k in 1:dim(ResultsMerged)[1]){
  Maximum_padj[k]  <- max(ResultsMerged$padj_TvsNT[k],ResultsMerged$padj_TvsH[k],ResultsMerged$padj_Cell[k]) 
}
#Replace NA for 1
Maximum_padj[which(is.na(Maximum_padj))] <- 1

#Create Color vector for 2d scatter plots
Scatter2DColor <- ifelse(Maximum_padj<0.05, rgb(150,150,150,150,maxColorValue=255), rgb(217,217,217,90,maxColorValue=255))
Scatter2DColor[which(ResultsMerged$TCGA_Cell_Target == "All_Upregulated")] <- rgb(20,20,20,255,maxColorValue=255)

########## FoldChange Correlations
x <- ResultsMerged$log2FoldChange_TvsNT
y <- ResultsMerged$log2FoldChange_TvsH

w <- round(cor(x,y),4)

png(file="./05_Results/Intersection/Correlations/Scatter_TvsNT-TvsH.png", res=350, height = 2000, width = 2000)
plot(x,y, xlim = c(-10,10), ylim=c(-10,10),
     main= "Correlation TvsNT-TvsH",
     xlab = "log2(Fold Change) - TvsNT",
     ylab = "log2(Fold Change) - TvsH",
     cex=0.8, cex.axis=1.4,
     cex.lab=1.5, cex.main=1.6,
     col=Scatter2DColor,
     #col=rgb(50,50,50,50,maxColorValue=255),
     pch=16)
abline(lm(y~x), 
       col="red", lwd=2, lty=4)

abline(h=1, lty=3, lwd=2, col=rgb(50,50,50,200,maxColorValue=255))
abline(v=1, lty=3, lwd=2, col=rgb(50,50,50,200,maxColorValue=255))

text(4.5, -8, paste("Pearson = ",w,sep=""), cex=1.2, font=3)
text(4.5, -10, paste("r-squared = ",
                      round(summary(lm(y~x))$adj.r.squared,4),
                      sep=""), cex=1.2, font=3)

dev.off()


########## FoldChange Correlations
x <- ResultsMerged$log2FoldChange_TvsNT
y <- ResultsMerged$log2FoldChange_Cell

w <- round(cor(x,y),4)

png(file="./05_Results/Intersection/Correlations/Scatter_TvsNT-Cell.png", res=350, height = 2000, width = 2000)
plot(x,y, xlim = c(-10,10), ylim=c(-10,10),
     main= "Correlation TvsNT-Cell",
     xlab = "log2(Fold Change) - TvsNT",
     ylab = "log2(Fold Change) - Cell",
     cex=0.8, cex.axis=1.4,
     cex.lab=1.5, cex.main=1.6,
     col=Scatter2DColor,
     pch=16)
abline(lm(y~x), 
       col="red", lwd=2, lty=4)
abline(h=1, lty=3, lwd=2, col=rgb(50,50,50,200,maxColorValue=255))
abline(v=1, lty=3, lwd=2, col=rgb(50,50,50,200,maxColorValue=255))

text(4.5, -8, paste("Pearson = ",w,sep=""), cex=1.2, font=3)
text(4.5, -10, paste("r-squared = ",
                     round(summary(lm(y~x))$adj.r.squared,4),
                     sep=""), cex=1.2, font=3)
dev.off()

########## FoldChange Correlations
x <- ResultsMerged$log2FoldChange_TvsH
y <- ResultsMerged$log2FoldChange_Cell

w <- round(cor(x,y),4)

png(file="./05_Results/Intersection/Correlations/Scatter_Cell-TvsH.png", res=350, height = 2000, width = 2000)
plot(x,y, xlim = c(-10,10), ylim=c(-10,10),
     main= "Correlation TvsH-Cell",
     xlab = "log2(Fold Change) - TvsH",
     ylab = "log2(Fold Change) - Cell",
     cex=0.8, cex.axis=1.4,
     cex.lab=1.5, cex.main=1.6,
     col=Scatter2DColor,
     pch=16)
abline(lm(y~x), 
       col="red", lwd=2, lty=4)

abline(h=1, lty=3, lwd=2, col=rgb(50,50,50,200,maxColorValue=255))
abline(v=1, lty=3, lwd=2, col=rgb(50,50,50,200,maxColorValue=255))

text(4.5, -8, paste("Pearson = ",w,sep=""), cex=1.2, font=3)
text(4.5, -10, paste("r-squared = ",
                     round(summary(lm(y~x))$adj.r.squared,4),
                     sep=""), cex=1.2, font=3)
dev.off()

# Creates the 3D scatterplot
# of log2FC correlation
# x, y and z coordinates
x <- ResultsMerged$log2FoldChange_TvsNT
y <- ResultsMerged$log2FoldChange_TvsH
z <- ResultsMerged$log2FoldChange_Cell

# Generate color status for graphs
ValueColor <- ifelse(Maximum_padj<0.05, 1, 0)
ValueColor[which(ResultsMerged$TCGA_Cell_Target == "All_Upregulated")] <- 2

#rgb(253,212,158,40,maxColorValue=255),
ColorSet <- c(
  rgb(35,139,69,30,maxColorValue=255),
  rgb(227,74,51,150,maxColorValue=255),
  rgb(122,1,119,150,maxColorValue=255)
)

ColorGrays <- c(
  rgb(217,217,217,30,maxColorValue=255),
  rgb(150,150,150,70,maxColorValue=255),
  rgb(49,49,49,255,maxColorValue=255)
)

# Special thanks to Alboukadel Kassambara and his blog
# "Statistical tools for high-throughput data analysis"
# For this part of the code
# http://www.sthda.com/english/wiki/impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization

#Function for scatter3d plotting with shadow projection
scatter3D_fancy <- function(x, y, z,..., colvar = ValueColor)
{
  panelfirst <- function(pmat) {
    
    XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = colvar, pch = ".", 
              cex = 2, add = TRUE, colkey = FALSE, col=ColorGrays)

    XY <- trans3D(x = rep(min(x), length(x)), y, z, pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = colvar, pch = ".", 
              cex = 2, add = TRUE, colkey = FALSE, col=ColorGrays)
  
    XY <- trans3D(x, y = rep(max(y), length(y)), z, pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = colvar, pch = ".", 
              cex = 2, add = TRUE, colkey = FALSE, col=ColorGrays)
    
    }
  scatter3D(x, y, z, ..., colvar = colvar, panel.first=panelfirst,
            colkey = FALSE) 
}

#Make the plot itself
png(paste("./05_Results/Intersection/Correlations/Scatter3d_projected",".png", sep=""), width=5000, height=7000, res=1200, type='cairo')
scatter3D_fancy(x, y, z, pch = 20, cex = 0.5,
                ticktype = "detailed", theta=60, phi = 45, d = 2,
                col=ColorSet,
                xlab="TvsNT",
                ylab="TvsH",
                zlab="Cell",
                ylim=c(-7.5,10))

dev.off()

####### TNBC Marker level on cell lines
MarkerList <- c("ESR1","ERBB2","PGR")

CellOrder <- c("HS578T_Daemen", 
               "HCC1806_Daemen",
               "HCC1937_Daemen", "HCC1937_Varley", 
               "HCC70_Daemen", "HCC70_Varley", 
               "HCC1143_Daemen", "HCC1143_Varley",
               "BT549_LNBio", "BT549_Daemen", "BT549_Varley", 
               "MDAMB157_Daemen", "MDAMB157_Varley",
               "MDAMB231_LNBio", "MDAMB231_Daemen", "MDAMB231_Varley", 
               "MDAMB436_LNBio", "MDAMB436_Varley",
               "MDAMB468_LNBio", "MDAMB468_Varley", 
               "MCF7_LNBio", "MCF7_Daemen", "MCF7_Varley",
               "ZR751_Daemen", "ZR751_Varley",
               "T47D_Daemen", "T47D_Varley",
               "SKBR3_LNBio", "SKBR3_Daemen", "SKBR3_Varley")

#
for(MarkerNow in MarkerList){
Data <- NormalizedGeneExpression[MarkerNow, CellOrder]
png(paste("./05_Results/CellLines/Marker_",MarkerNow,".png", sep=""), width=9000, height=3500, res=1200, type='cairo')
  par(mar=c(9,4,0.2,0.2), lwd = 2)
  barplot(as.numeric(Data),
        las=2, names=names(Data),
        col="white", lwd=2,
        ylab="NormRSEM")
dev.off()
}
