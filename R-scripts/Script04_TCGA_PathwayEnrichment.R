#
# Script04_TCGA_PathwayEnrichment.R
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
# This script intends to perform pathway enrichments using DE genes from TCGA cohort
#
# Inputs:
#		DE genes list
#
# Outputs:
#		Pathway enrichments
#

# Load libraries
library(methods)
library(RColorBrewer)
require(GO.db)
library(goseq)
library(clusterProfiler)
library(DOSE)
library(ReactomePA)
library(DO.db)
library(reactome.db)
library(VennDiagram)
library(org.Hs.eg.db)

#Create Results directory
dir.create("./04_Results", showWarnings = FALSE)
dir.create("./04_Results/GO", showWarnings = FALSE)
dir.create("./04_Results/REACTOME", showWarnings = FALSE)
dir.create("./04_Results/DO", showWarnings = FALSE)

#Drop SessionInfo about packages version
sink("./04_Results/sessionInfo.txt")
sessionInfo()
print("Date executed:")
print(Sys.time())
sink()

# Read previous DESeq2 output
TCGAExpression <- read.table(file="./03_Results/results_DESeq2_intersected_naming.txt", stringsAsFactors = FALSE)
TCGAExpression <- data.frame(TCGAExpression)

#################### GO ENRICHMENT
#################### START

# Get DE genes in both comparisons
TCGAExpression_DE <- TCGAExpression[which(TCGAExpression$BothTarget == "Both"), ]

# Create the Gene Vector
# Gene vector is a simple vector stating 1 for DE genes and
# 0 for non-de genes.
gene.vector_DE <- as.integer(TCGAExpression$geneSymbol_final%in%TCGAExpression_DE$geneSymbol_final)

# As this vector is unamed, put the right names on it
names(gene.vector_DE) <- TCGAExpression$geneSymbol_final

#Probability Weighting Function for Gene length
# Also plot the model
# the second position is for genome version and the third is for ID type
# hg19 is TCGA assembly version
# geneSymbol is intended to Gene Names
# knownGene is intended to GeneIDs
# To check all human genomes and supported ids, uncomment:
# supportedGenomes()[which(supportedGenomes()[,"species"] == "Human"),]
png(file=paste("./04_Results/GO/ProbabilityWeightingFunction",".png", sep=""), width=2000, height=2000, res = 400)
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
write.table(GO.wall_CC, file=paste("./04_Results/GO/CC",".txt",sep=""), sep="\t")
write.table(GO.wall_BP, file=paste("./04_Results/GO/BP",".txt",sep=""), sep="\t")
write.table(GO.wall_MF, file=paste("./04_Results/GO/MF",".txt",sep=""), sep="\t")

# Split GO table for each Go ontology and FDR values
GO.wall_fdr_CC <- GO.wall[which(GO.wall[,"ontology"]=='CC' & GO.wall[,"over_represented_fdr"] < 0.05),]
GO.wall_fdr_BP <- GO.wall[which(GO.wall[,"ontology"]=='BP' & GO.wall[,"over_represented_fdr"] < 0.05),]
GO.wall_fdr_MF <- GO.wall[which(GO.wall[,"ontology"]=='MF' & GO.wall[,"over_represented_fdr"] < 0.05),]

# Define category order
MainNames <- c("Cellular component", "Biological Process", "Molecular Funcion")

# Plot First 15 cat in each
i <- 0
for(pineapple in c("GO.wall_fdr_CC", "GO.wall_fdr_BP", "GO.wall_fdr_MF")){
  i <- i+1
  get(pineapple) -> PlotNowGO
  
  png(file=paste("./04_Results/GO/GOEnrichPlot_15first_",pineapple,".png",sep=""),
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
  
  png(file=paste("./04_Results/GO/GOEnrichPlot_all_",pineapple,".png",sep=""),
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
geneIDList_DE <- bitr(TCGAExpression_DE$geneSymbol_final, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb="org.Hs.eg.db")
geneIDList_Universe <- bitr(TCGAExpression$geneSymbol_final, fromType="SYMBOL",
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
write.table(REACTOME_Enriched_summary, file=paste("./04_Results/REACTOME/REAC.enrichment",".txt",sep=""), sep="\t")

#Get fdr <0.05
REACTOME_Enriched_fdr <- REACTOME_Enriched_summary[which(REACTOME_Enriched_summary$p.adjust<0.05),]

# REACTOME Plot
# First 15
REACTOME_Enriched_fdr[1:15,]-> PlotNowREACT

png(file=paste("./04_Results/REACTOME/REACT_EnrichPlot_15first.png",sep=""),
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

png(file=paste("./04_Results/REACTOME/REACT_EnrichPlot_All.png",sep=""),
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
png(file=paste("./04_Results/REACTOME/REACT_EnrichMap_NoNames.png",sep=""),
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
png(file=paste("./04_Results/REACTOME/REACT_EnrichMap.png",sep=""),
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
png(file=paste("./04_Results/REACTOME/REACT_cnetplot.png",sep=""),
    width=6000, heigh=6000, res=450)
cnetplot(REACTOME_Enriched, categorySize="Count", foldChange=NULL,
         showCategory = 1)
dev.off()

#################### REACTOME
#################### END

#################### Disease Ontology Semantic and Enrichment analysis
#################### START

# Another tool for conversion from gene name to gene ID, after HUGO translation
geneIDList_DE <- bitr(TCGAExpression_DE$geneSymbol_final, fromType="SYMBOL",
                      toType="ENTREZID", OrgDb="org.Hs.eg.db")
geneIDList_Universe <- bitr(TCGAExpression$geneSymbol_final, fromType="SYMBOL",
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
write.table(DO_Enriched_summary, file=paste("./04_Results/DO/DO.enrichment",".txt",sep=""), sep="\t")

#Get fdr <0.05
DO_Enriched_fdr <- DO_Enriched_summary[which(DO_Enriched_summary$p.adjust<0.05),]


# DO Plot
# First 15
DO_Enriched_fdr[1:15,]-> PlotNowDO

png(file=paste("./04_Results/DO/DO_EnrichPlot_15first.png",sep=""),
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

png(file=paste("./04_Results/DO/DO_EnrichPlot_All.png",sep=""),
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
