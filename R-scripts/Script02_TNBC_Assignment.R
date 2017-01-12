#
# Script02_TNBC_Assignment.R
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
# This script intends to assign breast marker status using RNA levels and select the best mixture model comparing to IHC data
#
# Inputs:
#		Non-normalized gene expression values (RSEM, RPKM, etc) from Script01_GDC-Download
#           and clinical data from GDC (nationwidechildrens.org_clinical_patient_brca.txt)
#
# Outputs:
#		Normalized gene expression values and classification table
#
# BRCA SDRF file could be obtained here:
# Clinical data could be obtained here:
# https://gdc-portal.nci.nih.gov/legacy-archive/files/735bc5ff-86d1-421a-8693-6e6f92055563
#

# Load libraries
library(methods)
library(mclust)
library(MineICA)
library(EBSeq)
library(RColorBrewer)

#Create Results directory
dir.create("./02_Results", showWarnings = FALSE)

#Drop SessionInfo about packages version
sink("./02_Results/sessionInfo.txt")
sessionInfo()
print("Date executed:")
print(Sys.time())
sink()

# Reads the count table
# The list must contain genes on rows and patients on cols
# Expression values should be raw (or Raw RSEMs, unormalized)
GeneExpression <- as.matrix(read.table("./01_Results/BRCA.PrimaryTumor.txt", stringsAsFactors = FALSE))
#Change dots for dashes again
colnames(GeneExpression) <- gsub("\\.","-", colnames(GeneExpression))

# Proccedes with Bullard upper quantile normalization using EBSeq function
# See 
#   Bullard, J. H., Purdom, E., Hansen, K. D., & Dudoit, S. (2010).
#   Evaluation of statistical methods for normalization
#   and differential expression in mRNA-Seq experiments.
#   BMC Bioinformatics, 11, 94. http://doi.org/10.1186/1471-2105-11-94

# It will generate an vector for normalization
BullardNormalization <- QuantileNorm(GeneExpression,.75)

# Get the normalization values and apply to the list
NormalizedGeneExpression <- GetNormalizedMat(GeneExpression, BullardNormalization)

# Write it down
write.table(NormalizedGeneExpression, file="./02_Results/BRCA.PrimaryTumor_normalized.txt")

# Normalization proccess checking
# Create an color ramp for all those lines
cols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(colnames(GeneExpression)))

#Plot Raw RSEM densities
png(file="./02_Results/RSEM_raw_plots_all.png", width=4000, height = 3000, res=600)
# Start the plot itself by the first column
plot(density(log2(GeneExpression[,1])),
     main= "Raw log2(RSEM) Density",
     ylab="Density",
     xlab="",
     col=cols[1],
     lwd=3,
     cex=1.5,
     cex.lab=1.5,
     yaxt="n",
     xaxt="n",
     bty="n",
     ylim=c(0,0.15)
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


#Plot Normalized RSEM densities
png(file="./02_Results/RSEM_Normalized_plots_all.png", width=4000, height = 3000, res=600)
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
     ylim=c(0,0.15)
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

# Sets the marker list by names on Gene Expression Data
MarkerList <- c("ERBB2|2064","ESR1|2099","PGR|5241")

# Remove the pipes for nomenclature
MarkerList -> MarkerName
for(i in 1:length(MarkerList)){
  strsplit(MarkerList, "[|]")[[i]][1] -> MarkerName[i]
}

# Perform some Previous plotting
i<-0
for(GeneNow in MarkerList){
i<-i+1
png(paste("./02_Results/Histogram_",MarkerName[i],"_Raw_log2.png",sep=""), width=2000, height=1500, res=300) 
hist(log2(GeneExpression[GeneNow,]+1), breaks = 50, main= paste(MarkerName[i]," Raw",sep=""), xlab="Expression log2(RSEM+1)", ylab="Density", cex.main=1.6, cex.lab=1.6, cex.axis=1.6)
dev.off()

png(paste("./02_Results/Histogram_",MarkerName[i],"_Normalized_log2.png",sep=""), width=2000, height=1500, res=300) 
hist(log2(NormalizedGeneExpression[GeneNow,]+1), breaks = 50, main= paste(MarkerName[i]," Normalized",sep=""), xlab="Expression log2(RSEM+1)", ylab="Density", cex.main=1.6, cex.lab=1.6, cex.axis=1.6)
dev.off()
}

# Apply Mclust individually
# Set the number of groups to be split and de model
# the model should be empirically determined
ModelList <- list(
  c("V","V","V"),
  c("E","V","V"),
  c("V","E","V"),
  c("V","V","E"),
  c("V","E","E"),
  c("E","V","E"),
  c("E","E","V"),
  c("E","E","E"))

# Create the vector for agreement results
ModelAgreementTriple <- rep(0, length(ModelList))

# Number of groups in each marker
NumberGroups <- rep(2, length(MarkerList))

#Loop for model evaluation
for(ModelNumber in 1:length(ModelList)){

  #Pick model name from the list
  ModelName <- ModelList[[ModelNumber]]

  # Create the dir for this model
  dir.create(paste("./02_Results/", paste(ModelName, collapse=""), sep=""), showWarnings = FALSE)
  
  # Result Table
  ResultTable <- matrix(ncol=(length(MarkerList)*4), nrow=length(colnames(GeneExpression)))

  # Put the colnames in rownames of Result Table
  rownames(ResultTable) <- colnames(NormalizedGeneExpression)
  colnames(ResultTable) <- 1:(length(MarkerList)*4)
  # Analysis itself
  # set up a outside counter
  i<-0
  k<-0
  for(GeneNow in MarkerList){
    i<-i+1
    #Perform the MClustModel
    MClustModel <- Mclust(log2(NormalizedGeneExpression[GeneNow,]+1),
                      G=NumberGroups[i],
                      modelName = ModelName[i])

    #Start plotting
    png(paste("./02_Results/", paste(ModelName, collapse=""),"/MClust_",MarkerName[i],"_%01d.png",sep=""), width=2000, height=1500, res=300)
    plotMix(mc=MClustModel,
        data=log2(NormalizedGeneExpression[GeneNow,]+1),
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
    GeneExpression[GeneNow,] -> ResultTable[,k]  
    colnames(ResultTable)[k] <- paste(MarkerName[i],"Raw",sep="_")
    k<-k+1
    NormalizedGeneExpression[GeneNow,] -> ResultTable[,k]  
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

    #Create the EBSeq Grouping scheme
    DESeq_Grouping <- matrix(ncol=1, nrow=length(rownames(ResultTable)))
    rownames(DESeq_Grouping) <- rownames(ResultTable)
    colnames(DESeq_Grouping) <- c("conditions")
    DESeq_Grouping[,1] <- ResultTable$TripleStatus

# Specific for Triple Negative Breast Cancer Markers
# Reads clinical file
clinical_patient <- read.table(file="./ExternalFiles/nationwidechildrens.org_clinical_patient_brca.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t", quote="\"")

# Remove lines one and two 
clinical_patient[c(-1,-2),] -> clinical_patient

#Add the column to the result table
ResultTable$er_status_by_ihc <- "NA"
ResultTable$pr_status_by_ihc <- "NA"
ResultTable$her2_status_by_ihc <- "NA"

# Loop for TCGA Barcode filling
for(w in 1:length(rownames(ResultTable))){
j <- 0
j <- match(paste(strsplit(rownames(ResultTable)[w],"")[[1]][1:12], collapse=""),
           clinical_patient$bcr_patient_barcode)
ResultTable$er_status_by_ihc[w] <- clinical_patient$er_status_by_ihc[j]
ResultTable$pr_status_by_ihc[w] <- clinical_patient$pr_status_by_ihc[j]
ResultTable$her2_status_by_ihc[w] <- clinical_patient$her2_status_by_ihc[j]
}

# Nomenclature arrangement
# We have a lot of "Not Evaluated", " Equivocal", "NA" in TCGA for all those markers
# Those lines just changes everything to "unknown".
gsub("[[]", "", ResultTable$er_status_by_ihc) -> ResultTable$er_status_by_ihc
gsub("[]]", "", ResultTable$er_status_by_ihc) -> ResultTable$er_status_by_ihc
gsub("Not Evaluated","Unknown",ResultTable$er_status_by_ihc) -> ResultTable$er_status_by_ihc
gsub("Indeterminate","Unknown",ResultTable$er_status_by_ihc) -> ResultTable$er_status_by_ihc
gsub("Equivocal","Unknown",ResultTable$er_status_by_ihc) -> ResultTable$er_status_by_ihc
gsub("Not Available","Unknown",ResultTable$er_status_by_ihc) -> ResultTable$er_status_by_ihc
ResultTable$er_status_by_ihc[is.na(ResultTable$er_status_by_ihc)] <- "Unknown"

gsub("[[]", "", ResultTable$pr_status_by_ihc) -> ResultTable$pr_status_by_ihc
gsub("[]]", "", ResultTable$pr_status_by_ihc) -> ResultTable$pr_status_by_ihc
gsub("Not Evaluated","Unknown",ResultTable$pr_status_by_ihc) -> ResultTable$pr_status_by_ihc
gsub("Indeterminate","Unknown",ResultTable$pr_status_by_ihc) -> ResultTable$pr_status_by_ihc
gsub("Equivocal","Unknown",ResultTable$pr_status_by_ihc) -> ResultTable$pr_status_by_ihc
gsub("Not Available","Unknown",ResultTable$pr_status_by_ihc) -> ResultTable$pr_status_by_ihc
ResultTable$pr_status_by_ihc[is.na(ResultTable$pr_status_by_ihc)] <- "Unknown"

gsub("[[]", "", ResultTable$her2_status_by_ihc) -> ResultTable$her2_status_by_ihc
gsub("[]]", "", ResultTable$her2_status_by_ihc) -> ResultTable$her2_status_by_ihc
gsub("Not Evaluated","Unknown",ResultTable$her2_status_by_ihc) -> ResultTable$her2_status_by_ihc
gsub("Indeterminate","Unknown",ResultTable$her2_status_by_ihc) -> ResultTable$her2_status_by_ihc
gsub("Equivocal","Unknown",ResultTable$her2_status_by_ihc) -> ResultTable$her2_status_by_ihc
gsub("Not Available","Unknown",ResultTable$her2_status_by_ihc) -> ResultTable$her2_status_by_ihc
ResultTable$her2_status_by_ihc[is.na(ResultTable$her2_status_by_ihc)] <- "Unknown"

# Changes 1 by Negative and 2 by Positive
gsub(1, "Negative",ResultTable$ESR1_classification) -> ResultTable$ESR1_classification
gsub(2, "Positive",ResultTable$ESR1_classification) -> ResultTable$ESR1_classification

gsub(1, "Negative",ResultTable$PGR_classification) -> ResultTable$PGR_classification
gsub(2, "Positive",ResultTable$PGR_classification) -> ResultTable$PGR_classification

gsub(1, "Negative",ResultTable$ERBB2_classification) -> ResultTable$ERBB2_classification
gsub(2, "Positive",ResultTable$ERBB2_classification) -> ResultTable$ERBB2_classification

# IHC Triple Status Assignment
ResultTable$triple_status_by_ihc <- NA

for(k in 1:length(ResultTable$triple_status_by_ihc)){
  if(ResultTable$er_status_by_ihc[k]=="Unknown"){
  } else {
    if(ResultTable$pr_status_by_ihc[k]=="Unknown"){
    } else {
      if(ResultTable$her2_status_by_ihc[k]=="Unknown"){
      } else {
        if(ResultTable$er_status_by_ihc[k]==ResultTable$pr_status_by_ihc[k] &&
        ResultTable$pr_status_by_ihc[k]==ResultTable$her2_status_by_ihc[k] &&
        ResultTable$her2_status_by_ihc[k]=="Negative"){
          ResultTable$triple_status_by_ihc[k] <- "Triple"
        } else {
          ResultTable$triple_status_by_ihc[k] <- "NonTriple"
        }
        }
    }
  }
}

# RNASeq - IHC Matching
# Estrogen Receptor
ResultTable$er_match <- NA
for(k in 1:length(ResultTable$er_status_by_ihc)){
  if(ResultTable$er_status_by_ihc[k]=="Unknown"){
  } else {
    ResultTable$er_match[k] <- ResultTable$er_status_by_ihc[k]==ResultTable$ESR1_classification[k]
  }
  }

# Progesterone Receptor
ResultTable$pr_match <- NA
for(k in 1:length(ResultTable$pr_status_by_ihc)){
  if(ResultTable$pr_status_by_ihc[k]=="Unknown"){
  } else {
    ResultTable$pr_match[k] <- ResultTable$pr_status_by_ihc[k]==ResultTable$PGR_classification[k]
  }
}

# Her2 Receptor
ResultTable$her2_match <- NA
for(k in 1:length(ResultTable$her2_status_by_ihc)){
  if(ResultTable$her2_status_by_ihc[k]=="Unknown"){
  } else {
    ResultTable$her2_match[k] <- ResultTable$her2_status_by_ihc[k]==ResultTable$ERBB2_classification[k]
  }
}

# Triple Status
ResultTable$triple_match <- NA
for(k in 1:length(ResultTable$triple_status_by_ihc)){
  if(isNA(ResultTable$triple_status_by_ihc[k])){
  } else {
    ResultTable$triple_match[k] <- ResultTable$triple_status_by_ihc[k]==ResultTable$TripleStatus[k]
  }
}

# Create barplot with concordances (in percentage) from RNAseq with IHC
ForBarPlot <- c(round(100*table(ResultTable$er_match)[2]/(table(ResultTable$er_match)[1]+table(ResultTable$er_match)[2]), 2),
round(100*table(ResultTable$pr_match)[2]/(table(ResultTable$pr_match)[1]+table(ResultTable$pr_match)[2]), 2),
round(100*table(ResultTable$her2_match)[2]/(table(ResultTable$her2_match)[1]+table(ResultTable$her2_match)[2]), 2),
round(100*table(ResultTable$triple_match)[2]/(table(ResultTable$triple_match)[1]+table(ResultTable$triple_match)[2]), 2))

# Barplot itself
# saves X positions for text
barplot(ForBarPlot)[,1] -> x_positions
dev.off()

# Plot open
png(width=1500, height = 2000, res=300, file=paste("./02_Results/", paste(ModelName, collapse=""), "/AgreementMarkers_Model_", paste(ModelName, collapse=""),".png",sep=""))
par(lwd = 3, mar=c(9,5,4,4))
barplot(ForBarPlot,
        names.arg=c("ER","PR","HER2","Triple"),
        las=2,
        col="#FECB92",
        border="#FDB462",
        main="Classification agreement",
        ylab="Concordance (%)",
        cex.names = 2,
        cex.axis= 2,
        lwd=2,
        cex.lab=2,
        ylim=c(0,110))
text(x_positions, ForBarPlot+6, labels = round(ForBarPlot,1), cex=2)
dev.off()

# Stacked bar plot - IHC
matrix(nrow=3, ncol=3) -> IHC.classification
rownames(IHC.classification) <- c("Negative", "Positive", "Unknown")
colnames(IHC.classification) <- c("ER","PR","HER2")
table(ResultTable$er_status_by_ihc) -> IHC.classification[,1]
table(ResultTable$pr_status_by_ihc) -> IHC.classification[,2]
table(ResultTable$her2_status_by_ihc) -> IHC.classification[,3]
IHC.classification <- IHC.classification[c(2,1,3),]

#Saves x positions for text
barplot(IHC.classification) -> x_positions
dev.off()

# Plot 1
png(width=1500, height = 2200, res=300, file=paste("./02_Results/", paste(ModelName, collapse=""), "/MarkerStatus_IHC_wNumb_", paste(ModelName, collapse=""),".png",sep=""))
par(lwd = 3, mar=c(10,6,5,5))
barplot(IHC.classification,
        las=2,
        col=c("#E2A4A4","#B2B0F7","#BFBFBF"),
        border=c("#CE6666","#7370F1","#999999"),
        main="Marker Status by IHC",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,1100))
text(x_positions, IHC.classification[1,]/2, labels = IHC.classification[1,], cex=2)
text(x_positions, IHC.classification[1,]+(IHC.classification[2,]/2), labels = IHC.classification[2,], cex=2)
text(x_positions, IHC.classification[1,]+IHC.classification[2,]+(IHC.classification[3,]/2), labels = IHC.classification[3,], cex=1.5)
dev.off()

# Plot 2
png(width=1500, height = 2200, res=300, file=paste("./02_Results/", paste(ModelName, collapse=""), "/MarkerStatus_IHC_woNumb_", paste(ModelName, collapse=""),".png",sep=""))
par(lwd = 3, mar=c(10,6,5,5))
barplot(IHC.classification,
        las=2,
        col=c("#E2A4A4","#B2B0F7","#BFBFBF"),
        border=c("#CE6666","#7370F1","#999999"),
        main="Marker Status by IHC",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,1100))
dev.off()

# Pie plot for histology
png(width=4000, height = 4000, res=400, file=paste("./02_Results/", paste(ModelName, collapse=""), "/PiePlot_IHC_", paste(ModelName, collapse=""),".png",sep=""))
par(lwd = 3, mar=c(10,10,10,10))
pie(table(ResultTable$triple_status_by_ihc, useNA = "always"),
    col=c("#F8766D","#00BFC4","#999999"),
    border=0,
    init.angle=180,
    edges=10000,
    cex=2,
    labels=paste(c("Non-Triple","Triple","Unknown"),table(ResultTable$triple_status_by_ihc, useNA = "always"),sep=" - "))
dev.off()

# Stacked bar plot - RNASeq
matrix(nrow=2, ncol=3) -> RNASeq.classification
rownames(RNASeq.classification) <- c("Negative", "Positive")
colnames(RNASeq.classification) <- c("ER","PR","HER2")
table(ResultTable$ESR1_classification) -> RNASeq.classification[,1]
table(ResultTable$PGR_classification) -> RNASeq.classification[,2]
table(ResultTable$ERBB2_classification) -> RNASeq.classification[,3]
RNASeq.classification <- RNASeq.classification[c(2,1),]

#Saves x positions for text
barplot(RNASeq.classification) -> x_positions
dev.off()

# Plot 1
png(width=1500, height = 2200, res=300, file=paste("./02_Results/", paste(ModelName, collapse=""), "/MarkerStatus_RNAseq_wNumb_", paste(ModelName, collapse=""),".png",sep=""))
par(lwd = 3, mar=c(10,6,5,5))
barplot(RNASeq.classification,
        las=2,
        col=c("#E2A4A4","#B2B0F7","#BFBFBF"),
        border=c("#CE6666","#7370F1","#999999"),
        main="Marker Status by RNASeq",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,1100))
text(x_positions, RNASeq.classification[1,]/2, labels = RNASeq.classification[1,], cex=2)
text(x_positions, RNASeq.classification[1,]+(RNASeq.classification[2,]/2), labels = RNASeq.classification[2,], cex=2)
dev.off()

# Plot 2
png(width=1500, height = 2200, res=300, file=paste("./02_Results/", paste(ModelName, collapse=""), "/MarkerStatus_RNAseq_woNumb_", paste(ModelName, collapse=""),".png",sep=""))
par(lwd = 3, mar=c(10,6,5,5))
barplot(RNASeq.classification,
        las=2,
        col=c("#E2A4A4","#B2B0F7","#BFBFBF"),
        border=c("#CE6666","#7370F1","#999999"),
        main="Marker Status by RNASeq",
        ylab="",
        cex.names = 2,
        cex.axis= 2.5,
        lwd=2,
        cex.lab=2,
        ylim=c(0,1100))
dev.off()

# Pie plot for RNAseq
png(width=4000, height = 4000, res=400, file=paste("./02_Results/", paste(ModelName, collapse=""), "/PiePlot_RNASeq_", paste(ModelName, collapse=""),".png",sep=""))
par(lwd = 3, mar=c(10,10,10,10))
pie(table(ResultTable$TripleStatus),
    col=c("#F8766D","#00BFC4","#999999"),
    border=0,
    init.angle=90,
    edges=10000,
    cex=2,
    labels=paste(c("Non-Triple","Triple"),table(ResultTable$TripleStatus),sep=" - "))
dev.off()

# Saves everything!
write.csv2(ResultTable, file=paste("./02_Results/", paste(ModelName, collapse=""), "/Marker_Result_", paste(ModelName, collapse=""),".csv",sep=""), quote = TRUE, eol = "\n", na = "NA")
write.table(ResultTable, file=paste("./02_Results/", paste(ModelName, collapse=""), "/Marker_Result_", paste(ModelName, collapse=""),".txt",sep=""), quote = TRUE, eol = "\n", na = "NA")
write.table(DESeq_Grouping, file=paste("./02_Results/", paste(ModelName, collapse=""), "/DESeq2_Groups", paste(ModelName, collapse=""),".txt",sep=""), quote = FALSE)

# Saves in a vector all the values
ModelAgreementTriple[ModelNumber] <- as.numeric(ForBarPlot[length(ForBarPlot)])
}


# Saves the Agreeement
names(ModelAgreementTriple) <- unlist(lapply(ModelList, function(w) (paste(w, collapse="")) ) )
sink(file="./02_Results/ModelAgreement.percentage.txt")
print(ModelAgreementTriple)
sink()

# Get the highest
BestModel <- names(ModelAgreementTriple)[match(max(ModelAgreementTriple), ModelAgreementTriple)]

# Copy the specific file to root, without the model name.
file.copy(paste("./02_Results/",BestModel,"/DESeq2_Groups",BestModel,".txt",sep=""),
          paste("./02_Results/","DESeq2_Groups.txt",sep=""))

# Orders the vector
BestModel <- BestModel[order(-BestModel)]

# Barplot itself
# saves X positions for text
barplot(BestModel)[,1] -> x_positions
dev.off()

# Plot open
png(width=2500, height = 2000, res=300, file=paste("./02_Results/", "BestModelAgreement_IHCvsRNASeq.png",sep=""))
par(lwd = 3, mar=c(9,5,4,4))
barplot(BestModel,
        names.arg=names(BestModel),
        las=2,
        col="#FECB92",
        border="#FDB462",
        main="Classification agreement",
        ylab="Concordance (%)",
        cex.names = 2,
        cex.axis= 2,
        lwd=2,
        cex.lab=2,
        ylim=c(0,110))
text(x_positions, BestModel+6, labels = round(BestModel,1), cex=2)
dev.off()
