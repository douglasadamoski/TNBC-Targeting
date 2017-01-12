#
# Script06_PubMed_Other_Graphs.R
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
# This script intends to perform access information about genes in PubMed, canSAR database druggability and plot some graphs from experimental data
#
# Inputs:
#		Gene lists and tables
#
# Outputs:
#		A lot of cool graphs
#

# Load libraries
library(RISmed)
library(RColorBrewer)
library(jsonlite)
library(tools)
library(httr)
library(calibrate)
library(ggplot2)
library(mclust)
library(MineICA)
library(methods)
library(survival)
library(survminer)
library(survMisc)
library(RTCGA.clinical)

# Create Results directory
dir.create("./06_Results", showWarnings = FALSE)
dir.create("./06_Results/PubMed", showWarnings = FALSE)
dir.create("./06_Results/DGIdb", showWarnings = FALSE)
dir.create("./06_Results/DGIdb/Partials", showWarnings = FALSE)
dir.create("./06_Results/canSAR", showWarnings = FALSE)
dir.create("./06_Results/Extras", showWarnings = FALSE)
dir.create("./06_Results/Clinical", showWarnings = FALSE)

#Drop SessionInfo about packages version
sink("./06_Results/sessionInfo.txt")
sessionInfo()
print("Date executed:")
print(Sys.time())
sink()

# Read TCGA categories
Categories <- data.frame(read.table(file="./02_Results/DESeq2_Groups.txt", stringsAsFactors = FALSE, sep=" ", header=TRUE))

# Read normalized table
NormalizedTCGAExpression <- read.table(file="./02_Results/BRCA.PrimaryTumor_normalized.txt", stringsAsFactors = FALSE)
colnames(NormalizedTCGAExpression) <- sapply(colnames(NormalizedTCGAExpression), function(w) {gsub("\\.","-",w)})

# Read qPCR GBP1 cell lines and plot
qPCR_GBP1 <- data.frame(read.table(file="./ExternalFiles/qPCR_GBP1_CellLines.txt", stringsAsFactors = FALSE, sep="\t"))

dodge <- position_dodge(width = 5)
limits <- aes(ymax = qPCR_GBP1$Upper,
              ymin = qPCR_GBP1$Lower)
ColorFill <- c(rep(rgb(0,186,56,150,maxColorValue=255), 2),
                rep(rgb(97,156,255,150,maxColorValue=255), 4))
ColorBorder <- c(rep(rgb(0,186,56,255,maxColorValue=255), 2),
               rep(rgb(97,156,255,255,maxColorValue=255), 4))
p <- ggplot(data = qPCR_GBP1,
            aes(x = factor(rownames(qPCR_GBP1), level=(rownames(qPCR_GBP1))),
                y = qPCR_GBP1$Mean))
png(filename="./06_Results/Extras/qPCR_GBP1_Cells.png", res=600, width=2000, height=2000)
p + geom_bar(stat = "identity", position = dodge,
             fill=ColorFill, colour=ColorBorder, size = 1) +
  geom_errorbar(limits, position = dodge, width = 0.25, size=0.9, colour="black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Read proliferation from cell lines and plot
Proliferation1 <- data.frame(read.table(file="./ExternalFiles/Proliferation_shRNA1.txt", stringsAsFactors = FALSE, sep="\t"))

limits <- aes(ymax = Proliferation1$Upper,
              ymin = Proliferation1$Lower)

ColorFill <- c(rep(rgb(97,156,255,150,maxColorValue=255), 8),
                 rep(rgb(0,186,56,150,maxColorValue=255), 3))
ColorBorder <- c(rep(rgb(97,156,255,255,maxColorValue=255), 8),
                 rep(rgb(0,186,56,255,maxColorValue=255), 3))
p <- ggplot(data = Proliferation1,
            aes(x = factor(rownames(Proliferation1), level=(rownames(Proliferation1))),
                y = Proliferation1$Mean))
png(filename="./06_Results/Extras/Proliferation1.png", res=600, width=3000, height=2000)
p + geom_bar(stat = "identity", position = dodge,
             fill=ColorFill, colour=ColorBorder, size = 1) +
  geom_errorbar(limits, position = dodge, width = 0.25, size=0.9, colour="black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# Read proliferation 2
Proliferation2 <- data.frame(read.table(file="./ExternalFiles/Proliferation_shRNA2.txt", stringsAsFactors = FALSE, sep="\t"))
limits <- aes(ymax = Proliferation2$Upper,
              ymin = Proliferation2$Lower)

ColorFill <- c(rep(rgb(97,156,255,150,maxColorValue=255), 8),
               rep(rgb(0,186,56,150,maxColorValue=255), 3))
ColorBorder <- c(rep(rgb(97,156,255,255,maxColorValue=255), 8),
                 rep(rgb(0,186,56,255,maxColorValue=255), 3))
p <- ggplot(data = Proliferation2,
            aes(x = factor(rownames(Proliferation2), level=(rownames(Proliferation2))),
                y = Proliferation2$Mean))
png(filename="./06_Results/Extras/Proliferation2.png", res=600, width=3000, height=2000)
p + geom_bar(stat = "identity", position = dodge,
             fill=ColorFill, colour=ColorBorder, size = 1) +
  geom_errorbar(limits, position = dodge, width = 0.25, size=0.9, colour="black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Single graph
JointProliferation <- data.frame(matrix(ncol=3, nrow=dim(Proliferation1)[1]+dim(Proliferation2)[1]))
colnames(JointProliferation) <- colnames(Proliferation1)

k <- 0
for(step in 1:dim(Proliferation1)[1]){
  k + 1 -> k
  rownames(JointProliferation)[k] <- paste(rownames(Proliferation1)[step],"-a",sep="")
  JointProliferation$Mean[k] <- Proliferation1$Mean[step]
  JointProliferation$Upper[k] <- Proliferation1$Upper[step]
  JointProliferation$Lower[k] <- Proliferation1$Lower[step]
    k + 1 -> k
  rownames(JointProliferation)[k] <- paste(rownames(Proliferation2)[step],"-b",sep="")
  JointProliferation$Mean[k] <- Proliferation2$Mean[step]
  JointProliferation$Upper[k] <- Proliferation2$Upper[step]
  JointProliferation$Lower[k] <- Proliferation2$Lower[step]
  }


limits <- aes(ymax = JointProliferation$Upper,
              ymin = JointProliferation$Lower)

ColorFill <- c(rep(rgb(97,156,255,150,maxColorValue=255), 16),
               rep(rgb(0,186,56,150,maxColorValue=255), 6))
ColorBorder <- c(rep(rgb(97,156,255,255,maxColorValue=255), 16),
                 rep(rgb(0,186,56,255,maxColorValue=255), 6))
p <- ggplot(data = JointProliferation,
            aes(x= sort(c(seq(from=1, to=dim(Proliferation1)[1]*2, by=2),
                 seq(from=1, to=dim(Proliferation1)[1]*2, by=2)+0.8)),
                #x = factor(rownames(JointProliferation), level=(rownames(JointProliferation))),
                y = JointProliferation$Mean, width = 0.65))
png(filename="./06_Results/Extras/ProliferationJoint.png", res=600, width=3000, height=2000)
p + geom_bar(stat = "identity", position = dodge,
             fill=ColorFill, colour=ColorBorder, size = 1) +
  geom_errorbar(limits, position = dodge, width = 0.2, size=0.75, colour="black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_continuous(breaks = c(seq(from=1, to=dim(Proliferation1)[1]*2, by=2)+0.4))
dev.off()


## Plot SEM not SD
replicates <- sqrt(3)
JointProliferation -> JointProliferationSEM
JointProliferationSEM$Upper <- JointProliferation$Mean+((JointProliferation$Upper - JointProliferation$Mean)/replicates)
JointProliferationSEM$Lower <-JointProliferation$Mean-((JointProliferation$Mean - JointProliferation$Lower)/replicates)

limits <- aes(ymax = JointProliferationSEM$Upper,
              ymin = JointProliferationSEM$Lower)

ColorFill <- c(rep(rgb(97,156,255,150,maxColorValue=255), 16),
               rep(rgb(0,186,56,150,maxColorValue=255), 6))
ColorBorder <- c(rep(rgb(97,156,255,255,maxColorValue=255), 16),
                 rep(rgb(0,186,56,255,maxColorValue=255), 6))

p <- ggplot(data = JointProliferationSEM,
            aes(x= sort(c(seq(from=1, to=dim(Proliferation1)[1]*2, by=2),
                          seq(from=1, to=dim(Proliferation1)[1]*2, by=2)+0.8)),
                #x = factor(rownames(JointProliferationSEM), level=(rownames(JointProliferationSEM))),
                y = JointProliferationSEM$Mean, width = 0.65))
png(filename="./06_Results/Extras/ProliferationJoint_SEM.png", res=600, width=3000, height=1750)
p + geom_bar(stat = "identity", position = dodge,
             fill=ColorFill, colour=ColorBorder, size = 1) +
  geom_errorbar(limits, position = dodge, width = 0.2, size=0.75, colour="black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_continuous(breaks = c(seq(from=1, to=dim(Proliferation1)[1]*2, by=2)+0.4))
dev.off()

# Read qPCR shRNA Efficiency
qPCR_shRNA_Eff <- data.frame(read.table(file="./ExternalFiles/qPCR_shRNA_Efficiency.txt", stringsAsFactors = FALSE, sep="\t"))

limits <- aes(ymax = qPCR_shRNA_Eff$Upper,
              ymin = qPCR_shRNA_Eff$Lower)
#Gera o arcabouço
p <- ggplot(data = qPCR_shRNA_Eff,
            aes(x = factor(rownames(qPCR_shRNA_Eff), level=(rownames(qPCR_shRNA_Eff))),
                y = qPCR_shRNA_Eff$Mean))
# Gira as labels
png(filename="./06_Results/Extras/qPCR_shRNA_Eff.png", res=600, width=1300, height=2000)
p + geom_bar(stat = "identity", position = dodge,
             fill="#bcbddc", colour="#756bb1", size = 1) +
  geom_errorbar(limits, position = dodge, width = 0.25, size=0.9, colour="black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Read qPCR EGF trial
qPCR_EGF_trial <- data.frame(read.table(file="./ExternalFiles/qPCR_EGF_trial.txt", stringsAsFactors = FALSE, sep="\t"))

limits <- aes(ymax = qPCR_EGF_trial$Upper,
              ymin = qPCR_EGF_trial$Lower)

ColorFill <- c(rep(rgb(0,186,56,150,maxColorValue=255), 4),
               rep(rgb(97,156,255,150,maxColorValue=255), 8))
ColorBorder <- c(rep(rgb(0,186,56,255,maxColorValue=255), 4),
                 rep(rgb(97,156,255,255,maxColorValue=255), 8))
p <- ggplot(data = qPCR_EGF_trial,
            aes(x = factor(rownames(qPCR_EGF_trial), level=(rownames(qPCR_EGF_trial))),
                y = qPCR_EGF_trial$Mean))
png(filename="./06_Results/Extras/qPCR_EGF_trial.png", res=600, width=3100, height=2000)
p + geom_bar(stat = "identity", position = dodge,
             fill=ColorFill, colour=ColorBorder, size = 1) +
  geom_errorbar(limits, position = dodge, width = 0.25, size=0.9, colour="black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

## Plot SEM not SD
replicates <- sqrt(3)
qPCR_EGF_trial -> qPCR_EGF_trialSEM
qPCR_EGF_trialSEM$Upper <- qPCR_EGF_trial$Mean+((qPCR_EGF_trialSEM$Upper - qPCR_EGF_trialSEM$Mean)/replicates)
qPCR_EGF_trialSEM$Lower <- qPCR_EGF_trial$Mean-((qPCR_EGF_trialSEM$Mean - qPCR_EGF_trialSEM$Lower)/replicates)


limits <- aes(ymax = qPCR_EGF_trialSEM$Upper,
              ymin = qPCR_EGF_trialSEM$Lower)

ColorFill <- c(rep(rgb(0,186,56,150,maxColorValue=255), 4),
               rep(rgb(97,156,255,150,maxColorValue=255), 8))
ColorBorder <- c(rep(rgb(0,186,56,255,maxColorValue=255), 4),
                 rep(rgb(97,156,255,255,maxColorValue=255), 8))
p <- ggplot(data = qPCR_EGF_trialSEM,
            aes(x = factor(rownames(qPCR_EGF_trialSEM), level=(rownames(qPCR_EGF_trialSEM))),
                y = qPCR_EGF_trialSEM$Mean))
png(filename="./06_Results/Extras/qPCR_EGF_trial_SEM.png", res=600, width=3100, height=2000)
p + geom_bar(stat = "identity", position = dodge,
             fill=ColorFill, colour=ColorBorder, size = 1) +
  geom_errorbar(limits, position = dodge, width = 0.25, size=0.9, colour="black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# Read qPCR GBP1 cell lines
qPCR_GBP1 <- data.frame(read.table(file="./ExternalFiles/qPCR_GBP1_CellLines.txt", stringsAsFactors = FALSE, sep="\t"))

dodge <- position_dodge(width = 5)
limits <- aes(ymax = qPCR_GBP1$Upper,
              ymin = qPCR_GBP1$Lower)
ColorFill <- c(rep(rgb(0,186,56,150,maxColorValue=255), 2),
               rep(rgb(97,156,255,150,maxColorValue=255), 4))
ColorBorder <- c(rep(rgb(0,186,56,255,maxColorValue=255), 2),
                 rep(rgb(97,156,255,255,maxColorValue=255), 4))
p <- ggplot(data = qPCR_GBP1,
            aes(x = factor(rownames(qPCR_GBP1), level=(rownames(qPCR_GBP1))),
                y = qPCR_GBP1$Mean))
png(filename="./06_Results/Extras/qPCR_GBP1_Cells.png", res=600, width=2000, height=2000)
p + geom_bar(stat = "identity", position = dodge,
             fill=ColorFill, colour=ColorBorder, size = 1) +
  geom_errorbar(limits, position = dodge, width = 0.25, size=0.9, colour="black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



# Read normalized table
NormalizedTCGAExpression <- read.table(file="./03_Results/BRCA.all.normalized.txt", stringsAsFactors = FALSE)

#
TripleStatusPatients <- read.table(file="./03_Results/BRCA.all.information.txt", stringsAsFactors = FALSE)
ColorVector <- gsub("NonTriple",rgb(0,186,56,150,maxColorValue=255),TripleStatusPatients$conditions)
ColorVector <- gsub("Triple",rgb(97,156,255,150,maxColorValue=255),ColorVector)
ColorVector <- gsub("Healthy",rgb(248,118,109,150,maxColorValue=255),ColorVector)

########## EGFR vs GBP1 Correlations
x <- as.numeric(log2(NormalizedTCGAExpression["EGFR|1956",]+1))
y <- as.numeric(log2(NormalizedTCGAExpression["GBP1|2633",]+1))

w <- round(cor(x,y),4)

png(file="./06_Results/Extras/Correlation_Patients_withHealthy.png", res=350, height = 2000, width = 2000)
plot(x,y, #xlim = c(-10,10), ylim=c(-10,10),
     main= "Correlation",
     xlab = "EGFR - log2(RSEM+1)",
     ylab = "GBP1 - log2(RSEM+1)",
     cex=0.8, cex.axis=1.4,
     cex.lab=1.5, cex.main=1.6,
     col=ColorVector,
     #col=rgb(50,50,50,50,maxColorValue=255),
     pch=16)
abline(lm(y~x), 
       col="red", lwd=2, lty=4)

text(15, 8, paste("Pearson = ",w,sep=""), cex=1.2, font=3)
text(15, 7, paste("r-squared = ",
                     round(summary(lm(y~x))$adj.r.squared,4),
                     sep=""), cex=1.2, font=3)
dev.off()

# Read normalized table
NormalizedTCGAExpression <- read.table(file="./02_Results/BRCA.PrimaryTumor_normalized.txt", stringsAsFactors = FALSE)

#
TripleStatusPatients <- read.table(file="./02_Results/DESeq2_Groups.txt", stringsAsFactors = FALSE)
ColorVector <- ifelse(TripleStatusPatients=="NonTriple",rgb(0,186,56,150,maxColorValue=255),rgb(97,156,255,150,maxColorValue=255))

########## EGFR vs GBP1 Correlations
x <- as.numeric(log2(NormalizedTCGAExpression["EGFR|1956",]+1))
y <- as.numeric(log2(NormalizedTCGAExpression["GBP1|2633",]+1))

w <- round(cor(x,y),4)

png(file="./06_Results/Extras/Correlation_Patients.png", res=350, height = 2000, width = 2000)
plot(x,y, #xlim = c(-10,10), ylim=c(-10,10),
     main= "Correlation",
     xlab = "EGFR - log2(RSEM+1)",
     ylab = "GBP1 - log2(RSEM+1)",
     cex=0.8, cex.axis=1.4,
     cex.lab=1.5, cex.main=1.6,
     col=ColorVector,
     #col=rgb(50,50,50,50,maxColorValue=255),
     pch=16)
abline(lm(y~x), 
       col="red", lwd=2, lty=4)

text(15, 8, paste("Pearson = ",w,sep=""), cex=1.2, font=3)
text(15, 7, paste("r-squared = ",
                  round(summary(lm(y~x))$adj.r.squared,4),
                  sep=""), cex=1.2, font=3)
dev.off()

# Read normalized table
NormalizedCellExpression <- read.table(file="./05_Results/CellLines/RSEM_CellLines.normalized.txt", stringsAsFactors = FALSE)
TripleStatusCells <- read.table(file="Classification_CellLines.txt", stringsAsFactors = FALSE)$conditions
ColorVector <- ifelse(TripleStatusCells=="NonTriple",rgb(0,186,56,150,maxColorValue=255),rgb(97,156,255,150,maxColorValue=255))



########## EGFR vs GBP1 Correlations
x <- as.numeric(log2(NormalizedCellExpression["EGFR",]+1))
y <- as.numeric(log2(NormalizedCellExpression["GBP1",]+1))

w <- round(cor(x,y),4)

png(file="./06_Results/Extras/Correlation_CellLines.png", res=350, height = 2000, width = 2000)
plot(x,y, #xlim = c(-10,10), ylim=c(-10,10),
     main= "Correlation",
     xlab = "EGFR - log2(RSEM+1)",
     ylab = "GBP1 - log2(RSEM+1)",
     cex=1.2, cex.axis=1.4,
     cex.lab=1.5, cex.main=1.6,
     #col=rgb(50,50,50,50,maxColorValue=255),
     col=ColorVector,
     pch=16)
abline(lm(y~x), 
       col="red", lwd=2, lty=4)

text(14, 2.5, paste("Pearson = ",w,sep=""), cex=1.2, font=3)
text(14, 1, paste("r-squared = ",
                  round(summary(lm(y~x))$adj.r.squared,4),
                  sep=""), cex=1.2, font=3)
dev.off()

########## EGFR vs GBP1 Correlations
x <- as.numeric(log2(NormalizedCellExpression["ELAVL1",]+1))
y <- as.numeric(log2(NormalizedCellExpression["GLS",]+1))

w <- round(cor(x,y),4)

png(file="./06_Results/Extras/ELAVL1vsGLS.png", res=350, height = 2000, width = 2000)
plot(x,y, xlim = c(5,18.2), #ylim=c(-10,10),
     main= "Correlation",
     xlab = "EGFR - log2(RSEM+1)",
     ylab = "GBP1 - log2(RSEM+1)",
     cex=1.2, cex.axis=1.4,
     cex.lab=1.5, cex.main=1.6,
     #col=rgb(50,50,50,50,maxColorValue=255),
     col=ColorVector,
     pch=16)
abline(lm(y~x), 
       col="red", lwd=2, lty=4)

textxy(x,y,colnames(NormalizedCellExpression), cex=0.4, offset=0.8)

text(14, 2.5, paste("Pearson = ",w,sep=""), cex=1.2, font=3)
text(14, 1, paste("r-squared = ",
                  round(summary(lm(y~x))$adj.r.squared,4),
                  sep=""), cex=1.2, font=3)
dev.off()

########## EGFR vs GBP1 Correlations
x <- as.numeric(log2(NormalizedCellExpression["EGFR",]+1))
y <- as.numeric(log2(NormalizedCellExpression["GBP1",]+1))

w <- round(cor(x,y),4)

png(file="./06_Results/Extras/Correlation_CellLines_names.png", res=350, height = 2000, width = 2000)
plot(x,y, xlim = c(5,18.2), #ylim=c(-10,10),
     main= "Correlation",
     xlab = "EGFR - log2(RSEM+1)",
     ylab = "GBP1 - log2(RSEM+1)",
     cex=1.2, cex.axis=1.4,
     cex.lab=1.5, cex.main=1.6,
     #col=rgb(50,50,50,50,maxColorValue=255),
     col=ColorVector,
     pch=16)
abline(lm(y~x), 
       col="red", lwd=2, lty=4)

textxy(x,y,colnames(NormalizedCellExpression), cex=0.4, offset=0.8)

text(14, 2.5, paste("Pearson = ",w,sep=""), cex=1.2, font=3)
text(14, 1, paste("r-squared = ",
                  round(summary(lm(y~x))$adj.r.squared,4),
                  sep=""), cex=1.2, font=3)
dev.off()


# Read intersection
ResultIntersected <- data.frame(read.table(file="./05_Results/Intersection/IntersectionAllTable.txt", stringsAsFactors = FALSE))

# Get up and down gene names
GeneListUp <- ResultIntersected$geneSymbol_final[ResultIntersected$TCGA_Cell_Target == "All_Upregulated"]
GeneListDown <- ResultIntersected$geneSymbol_final[ResultIntersected$TCGA_Cell_Target == "All_Downregulated"]

# Creates a table to get Pubmed results
resultsPubMed <- matrix(ncol=4, nrow=length(ResultIntersected$geneSymbol_final))
colnames(resultsPubMed) <- c("GeneAlone", "Cancer","Breast.Cancer", "TNBCancer")
rownames(resultsPubMed) <- ResultIntersected$geneSymbol_final

# Looping to grab paper number to every search string
# you should change "Sys.sleep(0.2)":
#
# ###
#
# IMPORTANT!!!!
#
# ###
#
# In order not to overload the E-utility servers,
# NCBI recommends that users post no more than three
# URL requests per second and limit large jobs to either
# weekends or between 9:00 PM and 5:00 AM
# Eastern time during weekdays. Failure to comply with
# this policy may result in an IP address being
# blocked from accessing NCBI.
#
# The Sys.sleep argument "stops" the script for some time (in sec)
# in order to reduce the amount of queries by second
#

DelayInSeconds <- 0.30

#External counter
i <- 0
for(geneNow in ResultIntersected$geneSymbol_final){
  i+1 -> i
  c(0,0,0,0) -> resultsPubMed[i,]
  #To bypass errors in 0 papers searches
  tryCatch({  
    
    # Perform pubmed query just for gene name
    res <- EUtilsSummary(paste(geneNow,sep=""),
                         type="esearch", db="pubmed", datetype='pdat',
                         mindate=1900, maxdate=2016, retmax=10000)
    QueryCount(res) -> resultsPubMed[i,1]
    Sys.sleep(DelayInSeconds)  
    
    # Perform pubmed query for gene name + cancer
    res <- EUtilsSummary(paste(geneNow," cancer",sep=""),
                         type="esearch", db="pubmed", datetype='pdat',
                         mindate=1900, maxdate=2016, retmax=10000)
    QueryCount(res) -> resultsPubMed[i,2]
    Sys.sleep(DelayInSeconds)
    
    # Perform pubmed query for gene name + breast cancer
    res <- EUtilsSummary(paste(geneNow," breast cancer",sep=""),
                         type="esearch", db="pubmed", datetype='pdat',
                         mindate=1900, maxdate=2016, retmax=10000)
    QueryCount(res) -> resultsPubMed[i,3]
    Sys.sleep(DelayInSeconds)
    
    # Perform pubmed query for gene name + triple negative breast cancer
    res <- EUtilsSummary(paste(geneNow," triple negative breast cancer",sep=""),
                         type="esearch", db="pubmed", datetype='pdat',
                         mindate=1900, maxdate=2016, retmax=10000)
    QueryCount(res) -> resultsPubMed[i,4]
    Sys.sleep(DelayInSeconds)
    
  }, error=function(e){})
  
}

# Write Result Table
write.table(resultsPubMed, file="./06_Results/PubMed/PubMedQuery_Table_AllGenes.txt")
resultsPubMed <- data.frame(read.table(file="./06_Results/PubMed/PubMedQuery_Table_AllGenes.txt", stringsAsFactors = FALSE))

###### PUBMED TARGETS
#Creates an vector for relation with TNBC
resultsPubMed$TNBC.Relationship <- 0

# Check if there are publications for triple negative breast cancer
# if you have at least one publication,change value to 1
for(position in 1:length(resultsPubMed$TNBC.Relationship)){
  if(resultsPubMed$TNBCancer[position]>=1){
    1 -> resultsPubMed$TNBC.Relationship[position]
  } else {
  }
}

#Creates an vector for relation with TNBC
resultsPubMed$BC.Relationship <- 0

# Check if there are publications for breast cancer
# if you have at least one publication,change value to 1
for(position in 1:length(resultsPubMed$BC.Relationship)){
  if(resultsPubMed$Breast.Cancer[position]>=1){
    1 -> resultsPubMed$BC.Relationship[position]
  } else {
  }
}

#Creates an vector for relation with TNBC
resultsPubMed$C.Relationship <- 0

# Check if there are publications for cancer
# if you have at least one publication,change value to 1
for(position in 1:length(resultsPubMed$C.Relationship)){
  if(resultsPubMed$Cancer[position]>=1){
    1 -> resultsPubMed$C.Relationship[position]
  } else {
  }
}

# Create a contingency table about the status
ContingencyTable <- matrix(nrow=9, ncol=2)
colnames(ContingencyTable) <- c("Unrelated", "Related")
rownames(ContingencyTable) <- c("AllGenes_TNBC","UpGenes_TNBC","DownGenes_TNBC",
                                "AllGenes_BreastCancer","UpGenes_BreastCancer","DownGenes_BreastCancer",
                                "AllGenes_Cancer","UpGenes_Cancer","DownGenes_Cancer")

data.frame(ContingencyTable) -> ContingencyTable

# Store tabulated results
table(resultsPubMed$TNBC.Relationship) -> ContingencyTable[1,]
table(resultsPubMed[GeneListUp,]$TNBC.Relationship) -> ContingencyTable[2,]
table(resultsPubMed[GeneListDown,]$TNBC.Relationship) -> ContingencyTable[3,]

table(resultsPubMed$BC.Relationship) -> ContingencyTable[4,]
table(resultsPubMed[GeneListUp,]$BC.Relationship) -> ContingencyTable[5,]
table(resultsPubMed[GeneListDown,]$BC.Relationship) -> ContingencyTable[6,]

table(resultsPubMed$C.Relationship) -> ContingencyTable[7,]
table(resultsPubMed[GeneListUp,]$C.Relationship) -> ContingencyTable[8,]
table(resultsPubMed[GeneListDown,]$C.Relationship) -> ContingencyTable[9,]

ContingencyTable$RelatedPercent <- 100*(ContingencyTable$Related/(ContingencyTable$Related+ContingencyTable$Unrelated))

# Write Result Table
write.table(ContingencyTable, file="./06_Results/PubMed/PubMedQuery_Contingency_Table.txt")

#BarPlot Genes
#Create the barplottable for stacked barplot output
BarPlotTable <- t(ContingencyTable[c(7,8,4,5,1,2),c(1,2)])

# Remove unrelated genes from subsequent categories
BarPlotTable[1,3] - BarPlotTable[1,1] -> BarPlotTable[1,3]
BarPlotTable[1,5] - BarPlotTable[1,1] - BarPlotTable[1,3] -> BarPlotTable[1,5]
BarPlotTable[1,4] - BarPlotTable[1,2] -> BarPlotTable[1,4]
BarPlotTable[1,6] - BarPlotTable[1,2] - BarPlotTable[1,4] -> BarPlotTable[1,6]

#Sum Unrelated and Related genes for each column
colSums(BarPlotTable) -> TotalGenesEach

#Stores the raw table
BarPlotTable -> BarPlotTable_RawNumbers

#Calculate the proportion for each column
BarPlotTable[1,]/TotalGenesEach -> BarPlotTable[1,]
BarPlotTable[2,]/TotalGenesEach -> BarPlotTable[2,]

#Transforms to percentage
BarPlotTable*100 -> BarPlotTable

#Get x positions of bars
x_positions <- barplot(BarPlotTable)
dev.off()

# Barplot with the numbers
png(width=2200, height = 2200, res=300, file="./06_Results/PubMed/PublicationEnrichment_up.png")
par(lwd = 3, mar=c(15,6,5,5))
barplot(BarPlotTable,
        las=2,
        col=c("lightgray", "#bcbddc"),
        border=c("black", "#756bb1"),
        #col=rev(brewer.pal(8,'Set3')[c(1,3)]),
        #c("#E2A4A4","#B2B0F7"),
        #border=c("#CE6666","#7370F1"),
        main="List enrichment for relevant genes",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,100))
text(x_positions, BarPlotTable[1,]/2, labels = BarPlotTable_RawNumbers[1,], cex=1.5)
text(x_positions, BarPlotTable[1,]+(BarPlotTable[2,]/2), labels = BarPlotTable_RawNumbers[2,], cex=1.5)

dev.off()


# Barplot without the numbers
png(width=2200, height = 2200, res=300, file="./06_Results/PubMed/PublicationEnrichment_wonumbers_up.png")
par(lwd = 3, mar=c(15,6,5,5))
barplot(BarPlotTable,
        las=2,
        col=c("lightgray", "#bcbddc"),
        border=c("black", "#756bb1"),
        #col=rev(brewer.pal(8,'Set3')[c(1,3)]),
        #c("#E2A4A4","#B2B0F7"),
        #border=c("#CE6666","#7370F1"),
        main="List enrichment for relevant genes",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,100))

dev.off()

#BarPlot Genes Downregulated
#Create the barplottable for stacked barplot output
BarPlotTable <- t(ContingencyTable[c(7,9,4,6,1,3),c(1,2)])


# Remove unrelated genes from subsequent categories
BarPlotTable[1,3] - BarPlotTable[1,1] -> BarPlotTable[1,3]
BarPlotTable[1,5] - BarPlotTable[1,1] - BarPlotTable[1,3] -> BarPlotTable[1,5]
BarPlotTable[1,4] - BarPlotTable[1,2] -> BarPlotTable[1,4]
BarPlotTable[1,6] - BarPlotTable[1,2] - BarPlotTable[1,4] -> BarPlotTable[1,6]

#Sum Unrelated and Related genes for each column
colSums(BarPlotTable) -> TotalGenesEach

#Stores the raw table
BarPlotTable -> BarPlotTable_RawNumbers

#Calculate the proportion for each column
BarPlotTable[1,]/TotalGenesEach -> BarPlotTable[1,]
BarPlotTable[2,]/TotalGenesEach -> BarPlotTable[2,]

#Transforms to percentage
BarPlotTable*100 -> BarPlotTable

#Get x positions of bars
x_positions <- barplot(BarPlotTable)
dev.off()

# Barplot with the numbers
png(width=2200, height = 2200, res=300, file="./06_Results/PubMed/PublicationEnrichment_down.png")
par(lwd = 3, mar=c(15,6,5,5))
barplot(BarPlotTable,
        las=2,
        col=c("lightgray", "#bcbddc"),
        border=c("black", "#756bb1"),
        #col=rev(brewer.pal(8,'Set3')[c(1,3)]),
        #c("#E2A4A4","#B2B0F7"),
        #border=c("#CE6666","#7370F1"),
        main="List enrichment for relevant genes",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,100))
text(x_positions, BarPlotTable[1,]/2, labels = BarPlotTable_RawNumbers[1,], cex=1.5)
text(x_positions, BarPlotTable[1,]+(BarPlotTable[2,]/2), labels = BarPlotTable_RawNumbers[2,], cex=1.5)

dev.off()


# Barplot without the numbers
png(width=2200, height = 2200, res=300, file="./06_Results/PubMed/PublicationEnrichment_wonumbers_down.png")
par(lwd = 3, mar=c(15,6,5,5))
barplot(BarPlotTable,
        las=2,
        col=c("lightgray", "#bcbddc"),
        border=c("black", "#756bb1"),
        #col=rev(brewer.pal(8,'Set3')[c(1,3)]),
        #c("#E2A4A4","#B2B0F7"),
        #border=c("#CE6666","#7370F1"),
        main="List enrichment for relevant genes",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,100))

dev.off()

# Pubmed barplot by class
BarPlotTable <- matrix(ncol=1, nrow=2)
rownames(BarPlotTable) <- c("ZeroOr1","2OrMore")

BarPlotTable[1,1] <- length(which(resultsPubMed[GeneListUp,"TNBCancer"] == 0 |
        resultsPubMed[GeneListUp,"TNBCancer"] == 1))

BarPlotTable[2,1] <- length(which(resultsPubMed[GeneListUp,"TNBCancer"] > 1))

BarPlotTable -> BarPlotTable_RawNumbers
BarPlotTable[1,1] <- 100*BarPlotTable_RawNumbers[1,1]/length(GeneListUp)
BarPlotTable[2,1] <- 100*BarPlotTable_RawNumbers[2,1]/length(GeneListUp)


# Barplot with the numbers
png(width=870, height = 2200, res=300, file="./06_Results/PubMed/PublicationEnrichment_Cat_wnumbers_down.png")
par(lwd = 3, mar=c(15,6,5,5))
barplot(BarPlotTable,
        las=2,
        col=c("lightgray", "#bcbddc"),
        border=c("black", "#756bb1"),
        #col=rev(brewer.pal(8,'Set3')[c(1,3)]),
        #c("#E2A4A4","#B2B0F7"),
        #border=c("#CE6666","#7370F1"),
        main="List enrichment for relevant genes",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,100))
text(x_positions, BarPlotTable[1,]/2, labels = BarPlotTable_RawNumbers[1,], cex=1.5)
text(x_positions, BarPlotTable[1,]+(BarPlotTable[2,]/2), labels = BarPlotTable_RawNumbers[2,], cex=1.5)

dev.off()


# Barplot without the numbers
png(width=870, height = 2200, res=300, file="./06_Results/PubMed/PublicationEnrichment_Cat_wonumbers_down.png")
par(lwd = 3, mar=c(15,6,5,5))
barplot(BarPlotTable,
        las=2,
        col=c("lightgray", "#bcbddc"),
        border=c("black", "#756bb1"),
        #col=rev(brewer.pal(8,'Set3')[c(1,3)]),
        #c("#E2A4A4","#B2B0F7"),
        #border=c("#CE6666","#7370F1"),
        main="List enrichment for relevant genes",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,100))

dev.off()

################################
# DGIdb query
# Get all possible interaction sources
writeBin(content(GET("http://dgidb.genome.wustl.edu/api/v1/interaction_sources.json",
      progress()), "raw"), paste("./06_Results/DGIdb/","DGIdbQuery.InteractionSources.json.txt",sep=""))
InteractionSources <- fromJSON(paste("./06_Results/DGIdb/","DGIdbQuery.InteractionSources.json.txt",sep=""))

#Get all possible drug types
writeBin(content(GET("http://dgidb.genome.wustl.edu/api/v1/drug_types.json",
      progress()), "raw"), paste("./06_Results/DGIdb/","DGIdbQuery.DrugTypes.json.txt",sep=""))
DrugTypes <- fromJSON(paste("./06_Results/DGIdb/DGIdbQuery.DrugTypes.json.txt",sep=""))

#Get all possible interaction types
writeBin(content(GET("http://dgidb.genome.wustl.edu/api/v1/interaction_types.json",
      progress()), "raw"), paste("./06_Results/DGIdb/","DGIdbQuery.InteractionTypes.json.txt",sep=""))
InteractionTypes <- fromJSON(paste("./06_Results/DGIdb/","DGIdbQuery.InteractionTypes.json.txt",sep=""))

#Remove spaces
InteractionTypes <- gsub(" ", "%20", InteractionTypes)

#Remove positive interactions
RemoveInteractions <- c ("activator", "agonist","inducer","potentiator","stimulator","positive%20allosteric%20modulator")
InteractionTypes <- InteractionTypes[!InteractionTypes %in% RemoveInteractions]

# Query the Database for our targets
writeBin(content(GET(
  paste("http://dgidb.genome.wustl.edu/api/v1/interactions.json?", 
        "genes=",
        paste(GeneListUp, collapse=","),
        "&drug_types=",
        paste(DrugTypes, collapse=","),
        "&interaction_sources=",
        paste(InteractionSources, collapse=","),
        "&interaction_types=",
        paste(InteractionTypes, collapse=","),
        sep=""),
progress()), "raw"), paste("./06_Results/DGIdb/","DGIdbQuery.Results.json.txt",sep=""))

# Read from JSON
DGIdbResult <- fromJSON(paste("./06_Results/DGIdb/","DGIdbQuery.Results.json.txt",sep=""), flatten=TRUE)

#Create the recipient table
DGIdb_summarized <- data.frame(matrix(ncol=0, nrow=length(DGIdbResult$matchedTerms$searchTerm)))

#Fill with some info
rownames(DGIdb_summarized) <- DGIdbResult$matchedTerms$searchTerm
DGIdb_summarized$geneName <- DGIdbResult$matchedTerms$geneName
DGIdb_summarized$geneLongName <- DGIdbResult$matchedTerms$geneLongName
DGIdb_summarized$amountOfDrugs <- "empty"
DGIdb_summarized$DrugList <- "empty"

# retrieve data from list
for(EachListPosition in 1:length(DGIdbResult$matchedTerms$interactions)){
  DGIdb_summarized$amountOfDrugs[EachListPosition] <- length(DGIdbResult$matchedTerms$interactions[[EachListPosition]]$drugName)
  DGIdb_summarized$DrugList[EachListPosition] <- paste(DGIdbResult$matchedTerms$interactions[[EachListPosition]]$drugName, collapse="|")
}

# Write table
write.table(DGIdb_summarized, file="./06_Results/DGIdb/DGIdb_summarized.txt")

# Loop for querying under DGIdb limit
# Query the Database for ALL targets (may take a longer time)
NumberOfRepeats <- ceiling(length(ResultIntersected$geneSymbol_final)/1000)

# Fill QueryRanges
QueryRanges <- matrix(ncol=2, nrow=NumberOfRepeats)
colnames(QueryRanges) <- c("start","end")
QueryRanges <- data.frame(QueryRanges)
QueryRanges$start <- seq(1, (NumberOfRepeats*1000), by=1000)
QueryRanges$end <- QueryRanges$start + 999

# Create the empty table to get the results
DGIdb_summarized_all <- DGIdb_summarized[NULL,]

# Loop for query
for(SliceNow in 1:NumberOfRepeats){
    writeBin(content(GET(
    paste("http://dgidb.genome.wustl.edu/api/v1/interactions.json?", 
          "genes=",
          paste(ResultIntersected$geneSymbol_final[QueryRanges$start[SliceNow]:QueryRanges$end[SliceNow]], collapse=","),
          "&drug_types=",
          paste(DrugTypes, collapse=","),
          "&interaction_sources=",
          paste(InteractionSources, collapse=","),
          "&interaction_types=",
          paste(InteractionTypes, collapse=","),
          sep=""),
    progress()), "raw"), paste("./06_Results/DGIdb/Partials/","DGIdbQuery.Results_all_Slice",SliceNow,".json.txt",sep=""))
  
  #Read from JSON
  DGIdbResult_temp <- fromJSON(paste("./06_Results/DGIdb/Partials/","DGIdbQuery.Results_all_Slice",SliceNow,".json.txt",sep=""), flatten=TRUE)
  
  #Create the recipient table
  DGIdb_summarized_temp <- data.frame(matrix(ncol=0, nrow=length(DGIdbResult_temp$matchedTerms$searchTerm)))
  
  #Fill with some info
  rownames(DGIdb_summarized_temp) <- DGIdbResult_temp$matchedTerms$searchTerm
  DGIdb_summarized_temp$geneName <- DGIdbResult_temp$matchedTerms$geneName
  DGIdb_summarized_temp$geneLongName <- DGIdbResult_temp$matchedTerms$geneLongName
  DGIdb_summarized_temp$amountOfDrugs <- "empty"
  DGIdb_summarized_temp$DrugList <- "empty"
  
  # retrieve data from list
  for(EachListPosition in 1:length(DGIdbResult_temp$matchedTerms$interactions)){
    tryCatch({ 
      DGIdb_summarized_temp$amountOfDrugs[EachListPosition] <- length(DGIdbResult_temp$matchedTerms$interactions[[EachListPosition]]$drugName)
      DGIdb_summarized_temp$DrugList[EachListPosition] <- paste(DGIdbResult_temp$matchedTerms$interactions[[EachListPosition]]$drugName, collapse="|")
      }, error=function(e){})
    }
 
  #Join table
  DGIdb_summarized_all <- rbind(DGIdb_summarized_all,DGIdb_summarized_temp)
  
  #Clear var
  remove(DGIdbResult_temp)
}

# Write table
write.table(DGIdb_summarized_all, file="./06_Results/DGIdb/DGIdb_summarized_all.txt")
data.frame(read.table(file="./06_Results/DGIdb/DGIdb_summarized_all.txt"))

# Summarize drugs
DrugSummary <- paste(DGIdb_summarized$DrugList, collapse="|")
DrugSummary <- unique(strsplit(DrugSummary, "\\|"))
DrugSummary_all <- paste(DGIdb_summarized_all$DrugList, collapse="|")
DrugSummary_all <- unique(strsplit(DrugSummary_all, "\\|"))

#Amount of Drugs
print(length(DrugSummary[[1]]))

# Make graph of targets with available drugs
BarPlotTable <- matrix(ncol=2, nrow=2)
colnames(BarPlotTable) <- c("All_Genes", "Targets")
rownames(BarPlotTable) <- c("DrugsUnavailable", "DrugsAvailable")

# Fill
BarPlotTable[1,1] <- length(ResultIntersected$geneSymbol_final)-dim(DGIdb_summarized_all)[1]
BarPlotTable[2,1] <- dim(DGIdb_summarized_all)[1]
BarPlotTable[1,2] <- length(GeneListUp)-dim(DGIdb_summarized)[1]
BarPlotTable[2,2] <- dim(DGIdb_summarized)[1]
BarPlotTable -> BarPlotTable_RawNumbers

#
BarPlotTable[1,] <- BarPlotTable[1,]/colSums(BarPlotTable_RawNumbers)
BarPlotTable[2,] <- BarPlotTable[2,]/colSums(BarPlotTable_RawNumbers)

#
BarPlotTable*100 -> BarPlotTable

# Barplot with the numbers
png(width=1300, height = 2200, res=300, file="./06_Results/DGIdb/DrugEnrichment_up.png")
par(lwd = 3, mar=c(15,6,5,5))
barplot(BarPlotTable,
        las=2,
        col=c("lightgray", "#bcbddc"),
        border=c("black", "#756bb1"),
        #col=rev(brewer.pal(8,'Set3')[c(1,3)]),
        #c("#E2A4A4","#B2B0F7"),
        #border=c("#CE6666","#7370F1"),
        main="List enrichment for relevant genes",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,100))
text(x_positions, BarPlotTable[1,]/2, labels = BarPlotTable_RawNumbers[1,], cex=1.5)
text(x_positions, BarPlotTable[1,]+(BarPlotTable[2,]/2), labels = BarPlotTable_RawNumbers[2,], cex=1.5)
dev.off()


##############
# canSAR
# As we can't find an API for canSAR, we grabbed the list
# manually from https://cansar.icr.ac.uk/cansar/tools/
# Annotation set. You should remove single and double quotes
canSAR_results <- read.table(file="./ExternalFiles/canSAR_output.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)

#
canSAR_results_subset <- canSAR_results[,c("search_term","number_of_3d_structure",
"number_of_3d_structure_druggable",
"ligand_druggability_score",
"ligand_druggability_percentile_rank")]

#
PubMedList <- GeneListUp[which(resultsPubMed[GeneListUp,"TNBCancer"] == 0 |
                                    resultsPubMed[GeneListUp,"TNBCancer"] == 1)]

#
rownames(canSAR_results_subset) <- canSAR_results_subset$search_term
canSAR_results_subset <- canSAR_results_subset[PubMedList,]

# Make graph of targets with available drugs
BarPlotTable <- c(length(GeneListUp),
dim(canSAR_results_subset)[1],
length(which(canSAR_results_subset$number_of_3d_structure > 0)),
length(which(canSAR_results_subset$number_of_3d_structure_druggable > 0)))
names(BarPlotTable) <- c("All","AvailableAt_canSAR", "WithStructure","WithDruggableStructure")

# Make more for plots
BarPlotTable -> BarPlotTable2
BarPlotTable -> BarPlotTable3

#Ordering
# Bar plot av/unav
BarPlotTable[1] <- BarPlotTable[1]-BarPlotTable[2]
BarPlotTable1 <- BarPlotTable[c(2,1)]

# Bar plot av/unav structure
BarPlotTable2[2] <- BarPlotTable2[2]-BarPlotTable2[3]
BarPlotTable2 <- BarPlotTable2[c(3,2)]

# Barplot drug non drug
BarPlotTable3[3] <- BarPlotTable3[3]-BarPlotTable3[4]
BarPlotTable3 <- BarPlotTable3[c(4,3)]

#
BarPlotTable <- matrix(ncol=3, nrow=2)
BarPlotTable[,1] <- BarPlotTable1
BarPlotTable[,2] <- BarPlotTable2
BarPlotTable[,3] <- BarPlotTable3

# Define x positions
x_positions <- barplot(BarPlotTable)
# Barplot with the numbers
png(width=1900, height = 2200, res=300, file="./06_Results/canSAR/canSAR_plot_wNumbers.png")
par(lwd = 3, mar=c(12,6,5,5))
barplot(BarPlotTable,
        las=2,
        col=c("#bcbddc","lightgray"),
        border=c("#756bb1","black"),
        #col=rev(brewer.pal(8,'Set3')[c(1,3)]),
        #c("#E2A4A4","#B2B0F7"),
        #border=c(,"#7370F1","#CE6666"),
        main="Druggability",
        names=c("All","WithStructure","WithDrugg"),
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2)
text(x_positions, BarPlotTable[1,]/2, labels = BarPlotTable[1,], cex=1.5)
text(x_positions, BarPlotTable[1,]+(BarPlotTable[2,]/2), labels = BarPlotTable[2,], cex=1.5)
dev.off()

# Save Table
write.csv2(canSAR_results_subset, file="./06_Results/canSAR/canSAR_results_subset.csv")
