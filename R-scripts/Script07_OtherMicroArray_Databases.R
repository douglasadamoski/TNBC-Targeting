#
# Script07_OtherMicroArray_Databases.R
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
# This script intends to access reprodutibiliy of TCGA breast cancer cohort in external microarray datasets
#
# Inputs:
#		ArrayExpress access codes list
#
# Outputs
#		A lot of TNBC classifications, DE genes tables and fancy graphs.
#

# Load libraries
library(limma)
library(ArrayExpress)
library(affy)
library(WGCNA)
library(annotate)
library(VennDiagram)
allowWGCNAThreads()
library(mclust)
library(MineICA)
library(methods)
library(RColorBrewer)

# If you need the agilent chip databases
# source("https://bioconductor.org/biocLite.R")
# biocLite(
# c("hgu133plus2.db","hgu133a.db","hgu133a2.db","hgu133a2cdf","hgu133a2frmavecs","hgu133a2probe","hgu133acdf","hgu133afrmavecs","hgu133aprobe","hgu133atagcdf","hgu133atagprobe","hgu133b.db","hgu133bcdf","hgu133bprobe","hgu133plus2.db","hgu133plus2cdf","hgu133plus2frmavecs","hgu133plus2probe","hgu219.db","hgu219cdf","hgu219probe")
# )

# Dir create
dir.create("./07_Results", showWarnings = FALSE)
dir.create("./07_Results/RawData", showWarnings = FALSE)
dir.create("./07_Results/Processed", showWarnings = FALSE)

#Drop SessionInfo about packages version
sink("./07_Results/sessionInfo.txt")
sessionInfo()
sink()

#Read RNAseq Results
RNASeqResults <- data.frame(read.table(file="./05_Results/Intersection/IntersectionAllTable.txt", stringsAsFactors = FALSE))

# Create dataset list from ArrayExpress
DatasetList <- c("E-MTAB-365","E-GEOD-65216","E-GEOD-12276","E-MTAB-1547","E-GEOD-4922","E-GEOD-3494","E-GEOD-1456")

# For loop to download the data and organize in folders
# will occupy HUGE space on disk! Take care!
for(DatasetNameNow in DatasetList){

  #Create raw data dir
  dir.create(paste("./07_Results/RawData/",DatasetNameNow,sep=""), showWarnings = FALSE)

  #Download data from ArrayExpress
  RawExperiment <- getAE(DatasetNameNow, type = "raw", path=paste("./07_Results/RawData/",DatasetNameNow,"/",sep=""))
}

#
# Special treatment for single datasets SDRFs
# Start

# E-MTAB-365 - Transcription profiling by array of
# breast cancer samples to define breast cancer subsets
DatasetList[1] -> DatasetNameNow
#Read the sample and data relationship format
SDRF <- read.delim(paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".sdrf.txt",sep=""),
                   check.names=FALSE,stringsAsFactors=FALSE)
#Select the CEL files from breast tissue
SDRF <- SDRF[SDRF[,"Characteristics [OrgansimPart]"] == "breast" & grepl("\\.CEL", SDRF[,"Array Data File"]),]
#
write.table(SDRF, paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".corrected.sdrf.txt",sep=""))

# E-GEOD-65216 - Expression profiling of breast
# cancer samples from Institut Curie (Maire cohort)
DatasetList[2] -> DatasetNameNow
#Read the sample and data relationship format
SDRF <- read.delim(paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".sdrf.txt",sep=""),
                   check.names=FALSE,stringsAsFactors=FALSE)
#Select the CEL files from breast tissue
SDRF <- SDRF[which(!SDRF$`FactorValue [sample_group]` == "CellLine" & !SDRF$`FactorValue [sample_group]` == "Healthy"),]
#
write.table(SDRF, paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".corrected.sdrf.txt",sep=""))
#

# E-GEOD-12276 - Transcription profiling by array
# of human primary breast tumor tissue
DatasetList[3] -> DatasetNameNow
#Read the sample and data relationship format
SDRF <- read.delim(paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".sdrf.txt",sep=""),
                   check.names=FALSE,stringsAsFactors=FALSE)
# Nothing to be selected
write.table(SDRF, paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".corrected.sdrf.txt",sep=""))

# E-MTAB-1547 - Translation profiling of human
# inflammatory breast cancer samples to investigate
# gene expression
DatasetList[4] -> DatasetNameNow
#Read the sample and data relationship format
SDRF <- read.delim(paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".sdrf.txt",sep=""),
                   check.names=FALSE,stringsAsFactors=FALSE)
#Nothing to be done
write.table(SDRF, paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".corrected.sdrf.txt",sep=""))

# E-GEOD-4922 - Transcription profiling of human breast
# cancer tumor samples from Uppsala and Singapore cohorts
DatasetList[5] -> DatasetNameNow
#Read the sample and data relationship format
SDRF <- read.delim(paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".sdrf.txt",sep=""),
                   check.names=FALSE,stringsAsFactors=FALSE)
#Select the CEL files from breast tissue and HG-U133A
SDRF <- SDRF[which(SDRF$`Characteristics [OrganismPart]` == "mammary gland" & SDRF$`Comment [Array Design URI]` =="http://www.ebi.ac.uk/aerep/result?queryFor=PhysicalArrayDesign&aAccession=A-AFFY-33"),]
#
write.table(SDRF, paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".corrected.sdrf.txt",sep=""))
#

# E-GEOD-3494 - Transcription profiling of human breast cancer
# samples to generate an expression signature for p53 status
DatasetList[6] -> DatasetNameNow
#Read the sample and data relationship format
SDRF <- read.delim(paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".sdrf.txt",sep=""),
                   check.names=FALSE,stringsAsFactors=FALSE)
#Select the CEL files from HG-U133A
SDRF <- SDRF[which(SDRF$`Comment [Array Design URI]` =="http://www.ebi.ac.uk/aerep/result?queryFor=PhysicalArrayDesign&aAccession=A-AFFY-33"),]

#
write.table(SDRF, paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".corrected.sdrf.txt",sep=""))

# E-GEOD-1456 - Transcription profiling of breast cancer tissue in a
# large population-based cohort of Swedish patients
DatasetList[7] -> DatasetNameNow
#Read the sample and data relationship format
SDRF <- read.delim(paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".sdrf.txt",sep=""),
                   check.names=FALSE,stringsAsFactors=FALSE)
#Select the CEL files from HG-U133A
SDRF <- SDRF[which(SDRF$`Comment [Array Design URI]` =="http://www.ebi.ac.uk/aerep/result?queryFor=PhysicalArrayDesign&aAccession=A-AFFY-33"),]
#
write.table(SDRF, paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".corrected.sdrf.txt",sep=""))
#

# End
# Special treatment for single datasets SDRFs
#

# For loop to analyze the data
# will take HUGE time and memory! Take care!
for(DatasetNameNow in DatasetList){

#Read the sample and data relationship format already corrected
SDRF <- read.table(paste("./07_Results/RawData/",DatasetNameNow,"/",DatasetNameNow,".corrected.sdrf.txt",sep=""),
                   check.names=FALSE,stringsAsFactors=FALSE)

# Read CEL files
x <- justRMA(filenames=paste("./07_Results/RawData/",DatasetNameNow,"/",SDRF[,"Array Data File"],sep=""),
                normalize = TRUE, background = TRUE, verbose=TRUE, destructive = TRUE)

# Get the affymetrix platform name
x@annotation -> platformName
paste(platformName,".db", sep="") -> platformName_pkg
library(platformName_pkg, character.only = TRUE)

# Look up the Gene Symbol, name, and Ensembl Gene ID for each of those IDs
ID <- featureNames(x)
Symbol <- getSYMBOL(ID, platformName_pkg)
Name <- as.character(lookUp(ID, platformName_pkg, "GENENAME"))

#Background correct and normalize
w <- backgroundCorrect(x,method="normexp")
y <- normalizeBetweenArrays(w,method="quantile")

# WGCNA collapseRows
y_genes <- collapseRows(y, rowGroup=Symbol, rowID=featureNames(x), method="MaxMean")

#Grab the expression data
GeneExpression <- data.frame(y_genes$datETcollapsed)

#Create processed folder
dir.create(paste("./07_Results/Processed/",DatasetNameNow,sep=""), showWarnings = FALSE)

#Create an color ramp for all those lines
cols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(colnames(GeneExpression)))
#Plot Raw RSEM densities
#Open PNG file. MAY BE A TROUBLE in RScripts running on linux sessions without X
#One may find usefull change it for PDF format
png(file=paste("./07_Results/Processed/",DatasetNameNow,"/NormIntensityDensity.png",sep=""), width=4000, height = 3000, res=600)
# Start the plot itself by the first column
plot(density(log2(GeneExpression[,1]+1)),
     main= "Normalized Intensity Density",
     ylab="Density",
     xlab="MA Intensity Density (log2)",
     col=cols[1],
     lwd=3,
     cex=1.5,
     cex.lab=1.5,
     yaxt="n",
     xaxt="n",
     bty="n",
     ylim=c(0,(max(density(log2(GeneExpression[,1]+1))$y)*1.1))
)

#Add all other columns, one by one
for(ActualPatient in 2:length(colnames(GeneExpression))){
  lines(density(log2(GeneExpression[,ActualPatient]+1)),
        col=cols[ActualPatient],
        lwd=3)  
}

# Add the axis
axis(side = 1, lwd = 2,cex.axis=1.4)
axis(side = 2, lwd = 2,cex.axis=1.4)
box(lwd=2)
dev.off()

# Sets the marker list by names on Gene Expression Data
MarkerName <- c("ERBB2","ESR1","PGR")

# Perform some Previous plotting
i<-0
for(GeneNow in MarkerName){
  i<-i+1
  png(paste("./07_Results/Processed/",DatasetNameNow,"/Histogram_",MarkerName[i],"_log2Intensity.png",sep=""), width=2000, height=1500, res=300) 
  hist(as.numeric(log2(GeneExpression[GeneNow,]+1)), breaks = 50, main= paste(MarkerName[i]," Normalized",sep=""), xlab="Expression log2(Intensity)", ylab="Density", cex.main=1.6, cex.lab=1.6, cex.axis=1.6)
  dev.off()
  }

###########
# mclust START
# the model was already selected and tested over TCGA data
ModelName <- c("E","E","E")

# Number of groups in each marker
NumberGroups <- rep(2, length(MarkerName))

# Result Table
ResultTable <- matrix(ncol=(length(MarkerName)*4), nrow=length(colnames(GeneExpression)))

# Put the colnames in rownames of Result Table
rownames(ResultTable) <- colnames(GeneExpression)
colnames(ResultTable) <- 1:(length(MarkerName)*4)
# Analysis itself
# set up a outside counter
i<-0
k<-0
for(GeneNow in MarkerName){
  i<-i+1
  #Perform the MClustModel
  MClustModel <- Mclust(as.numeric(log2(GeneExpression[GeneNow,]+1)),
                        G=NumberGroups[i],
                        modelName = ModelName[i])
  
  #Start plotting
  png(paste("./07_Results/Processed/",DatasetNameNow,"/MClust_",MarkerName[i],"_%01d.png",sep=""), width=2000, height=1500, res=300)
  plotMix(mc=MClustModel,
          data=as.numeric(log2(GeneExpression[GeneNow,]+1)),
          nbBreaks = 50,
          traceDensity = TRUE,
          title =MarkerName[i],
          cex.main=1.6,
          cex.lab=1.6,
          cex.axis=1.6,
          lwd=2,
          cex=1)
  
  plot(MClustModel, what = "BIC")
  plot(MClustModel, what = "classification")
  plot(MClustModel, what = "uncertainty")
  plot(MClustModel, what = "density")
  dev.off()
  
  # Save values in our premade table
  k<-k+1
  as.numeric(GeneExpression[GeneNow,]) -> ResultTable[,k]  
  colnames(ResultTable)[k] <- paste(MarkerName[i],"Raw",sep="_")
  k<-k+1
  as.numeric(GeneExpression[GeneNow,]) -> ResultTable[,k]  
  colnames(ResultTable)[k] <- paste(MarkerName[i],"Normalized",sep="_")  
  k<-k+1
  MClustModel$classification -> ResultTable[,k]
  colnames(ResultTable)[k] <- paste(MarkerName[i],"classification",sep="_") 
  k<-k+1
  as.numeric(MClustModel$uncertainty) -> ResultTable[,k]
  colnames(ResultTable)[k] <- paste(MarkerName[i],"uncertainty",sep="_")  
}

#Move to a dataframe
data.frame(ResultTable) -> ResultTable

#If you are using just breast cancer common markers, asssign the triple ones
ResultTable$TripleStatus <- "Missing"
for(ActualPatient in 1:length(ResultTable[,1])){
  
  if(rowSums(ResultTable[ActualPatient,c(3,7,11)]) == 3){
    ResultTable$TripleStatus[ActualPatient] <- "Triple"
  } else {
    ResultTable$TripleStatus[ActualPatient] <- "NonTriple"
  }
}

#Write table results
write.table(ResultTable, paste("./07_Results/Processed/",DatasetNameNow,"/mclustResult_Table.txt",sep=""))
write.table(GeneExpression, paste("./07_Results/Processed/",DatasetNameNow,"/GeneExpression_Table.txt",sep=""))

# Define table for barplot
MAClassification <- matrix(ncol=4,nrow=2)
MAClassification[,1] <- table(ResultTable$ESR1_classification)
MAClassification[,2] <- table(ResultTable$PGR_classification)
MAClassification[,3] <- table(ResultTable$ERBB2_classification)
MAClassification[,4] <- table(ResultTable$TripleStatus)
colnames(MAClassification) <- c("ESR1","PGR","ERBB2","Triple")

# Plot 1
x_positions <- barplot(MAClassification)
dev.off()
png(width=1600, height = 2200, res=300, file=paste("./07_Results/Processed/",DatasetNameNow,"/ClassificationPlot.png",sep=""))
par(lwd = 3, mar=c(10,6,5,5))
barplot(MAClassification,
        las=2,
        col=c("#E2A4A4","#B2B0F7","#BFBFBF"),
        border=c("#CE6666","#7370F1","#999999"),
        main=paste("Marker status - ",DatasetNameNow,sep=""),
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2)
text(x_positions, MAClassification[1,]/2, labels = MAClassification[1,], cex=2)
text(x_positions, MAClassification[1,]+(MAClassification[2,]/2), labels = MAClassification[2,], cex=2)
dev.off()

# Define design table
design <- model.matrix(~-1 + factor(ifelse(ResultTable$TripleStatus == "NonTriple",1,2)))
colnames(design) <- c("NonTriple","Triple")

# Fit
fit <- lmFit(GeneExpression,design)
fit <- eBayes(fit,trend=TRUE)
# Export table
completeTable <- topTable(fit,coef=2,number=dim(fit$coefficients)[1])
completeTable$Status <- "Irrelevant"
completeTable$Status[which(completeTable$logFC > 1 & completeTable$adj.P.Val < 0.001)] <- "Target"

#write table
write.table(completeTable,file=paste("./07_Results/Processed/",DatasetNameNow,"/limmaDEGenesResults.txt",sep=""))

# Define intersection between results
MA_DEGenes <- rownames(completeTable)[which(completeTable$Status == "Target")]
RNAseq_DEGenes <- RNASeqResults$geneSymbol_final[which(RNASeqResults$TCGA_Cell_Target == "All_Upregulated")]
UpregulatedGenes <- list(RNAseq=RNAseq_DEGenes,Marray=MA_DEGenes)

# Venn plot
venn.plot <- venn.diagram(UpregulatedGenes,
  height = 3000, width = 3000, resolution = 500, imagetype = "png",
  main=DatasetNameNow,
  main.cex=1,
  print.mode="raw",
  lwd=c("4","4"),
  col=c("#F8766D"),
  fill=c("#F8766D"),
  alpha=0.5,
  cat.cex = 1.5,
  cex=c(2,2,2),
  margin = 0.12,
  cat.dist=c(0.06,0.06),
  category.names = c("RNAseq","Marray"),
  filename=paste("./07_Results/Processed/",DatasetNameNow,"/VennDiagram_up.png",sep=""))

#Boxplotting DE Genes
#Create processed folder
dir.create(paste("./07_Results/Processed/",DatasetNameNow,"/Boxplots",sep=""), showWarnings = FALSE)

# Set the levels
mylevels <- c("NonTriple", "Triple")

# Create table just with DE and with FC cutoff
resultadosDE <- GeneExpression[RNAseq_DEGenes,]

#Define targets
ActualGeneList <- rownames(resultadosDE)
ActualGeneList <- ActualGeneList[which(!grepl("NA\\.",ActualGeneList))]
ActualGeneList <- ActualGeneList[which(!ActualGeneList == "NA")]

#Empties k value
k <- 0
# Boxplot Loop START >>>>>>>>>>>>>>>>
for(Target in ActualGeneList){
  k <- k+1

  # Stores the value in our temporary table
  PlottingNow <- data.frame(matrix(nrow=dim(GeneExpression)[2], ncol=0))
  rownames(PlottingNow) <- colnames(GeneExpression)
  PlottingNow$value <- as.numeric(log2(GeneExpression[Target,]+1))
  PlottingNow$condition <- ResultTable$TripleStatus
  
  # Perform the boxplot itself
  png(paste("./07_Results/Processed/",DatasetNameNow,"/Boxplots/",Target,"_",DatasetNameNow,".png",sep=""), width=5000, height=7000, res=1200, type='cairo')
  boxplot(PlottingNow$value~PlottingNow$condition,
          las = 1,
          col = c("#80D99A","#79D8DB"),
          border = c("#3AC664","#00BFC4"),
          outline=TRUE,
          range=0.5,
          boxwex=0.8,
          notch = TRUE,
          lwd=2,
          #width=levelProportions,
          ylab ="log2(Intensity+1)",
          main=Target,
          par(cex.lab=1.5, cex.axis=1),
          outpch=NA)
  
  #Adiciona os pontos posteriormente pra ficar show
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

#Finishes everything
}
