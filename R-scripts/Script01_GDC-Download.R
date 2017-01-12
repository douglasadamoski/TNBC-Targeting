#
# Script01_GDC-Download.R
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
# This script intends to download TCGA files from GDC data portal
#
# Inputs:
#		Manifest from GDC Legacy Data Portal and SDRF file
#
# Outputs:
#		Donwloaded Raw Gene Expression Values for legacy TCGA Gene Level
#
# BRCA SDRF file could be obtained here:
#   https://gdc-portal.nci.nih.gov/legacy-archive/files/eaa28346-539e-499a-8fa9-562f0401bd8b
#

# Load libraries
library(methods)
library(tools)
library(httr)

# Create 01_Results directory
dir.create("./01_Results", showWarnings = FALSE)
dir.create("./01_Results/01_Download", showWarnings = FALSE)

# Creates the session info
sink(file="./01_Results/sessionInfo.txt")
sessionInfo()
print("Date executed:")
print(Sys.time())
sink()

# Read the SDRF
SDRF <- data.frame(read.table(file="./ExternalFiles/unc.edu_BRCA.IlluminaHiSeq_RNASeqV2.1.12.0.sdrf.txt", stringsAsFactors = FALSE, header=TRUE, sep="\t"))

#Read the MANIFEST
# For our case we build the MANIFEST going to
# https://gdc-portal.nci.nih.gov/legacy-archive
# Under cases tab:
# 	Primary Site: Breast
# 	Cancer Program: TCGA
# 	Project: TCGA-BRCA
# 	Disease Type: Breast Invasive Carcinoma
# Under files tab:
#	  Data Category: Gene expression
#	  Data Type: Gene expression quantification
#	  Experimental Strategy: RNA-Seq
#	  Data Format: TXT
#	  Platform: Illumina HiSeq
#	  Access Level: open
# And download the MANIFEST file. Change the name here if needed
MANIFEST <- data.frame(read.table(file="./ExternalFiles/gdc_manifest.2016-09-27T12-16-37.965495.tsv",
                                  stringsAsFactors = FALSE, header=TRUE, sep="\t"))

#Select lines with ".rsem.genes.results" files
MANIFEST <- MANIFEST[which(grepl(".rsem.genes.results", MANIFEST$filename)),]

#Loop for file download and md5sum check
for(variable in 1:dim(MANIFEST)[1]){
  # Downloads the file from server
  print(MANIFEST$file[variable])
  writeBin(content(GET(
    paste("https://gdc-api.nci.nih.gov/legacy/data/", MANIFEST$id[variable], sep=""), progress()),
    "raw"),paste("./01_Results/01_Download/",MANIFEST$file[variable],sep=""))
  # Check the md5sum. If fails, download again.
  print(paste("MD5 Check:", (md5sum(paste("./01_Results/01_Download/",MANIFEST$file[variable],sep="")) ==  MANIFEST$md5[variable])))
  
  #Keep doing until get same md5sum
  while ( (md5sum(paste("./01_Results/01_Download/",MANIFEST$file[variable],sep="")) ==  MANIFEST$md5[variable]) == FALSE){
    #Caches the file from server
    writeBin(content(GET(
      paste("https://gdc-api.nci.nih.gov/legacy/data/", MANIFEST$id[variable], sep=""), progress()),
      "raw"),paste("./01_Results/01_Download/",MANIFEST$file[variable],sep=""))
    print(paste("MD5 Check:", (md5sum(paste("./01_Results/01_Download/",MANIFEST$file[variable],sep="")) ==  MANIFEST$md5[variable])))
    
  }
}

#Starts the table assembly by first file
Base <- read.table(paste("./01_Results/01_Download/",MANIFEST$file[1],sep=""), header=TRUE, stringsAsFactors = FALSE, sep="\t")
rownames(Base) <- Base[,1]
Base <- data.frame(Base)
Base <- Base[,NULL]

#Loop for table assembly
for(variable in 1:dim(MANIFEST)[1]){
  Base[,MANIFEST$file[variable]] <- read.table(paste("./01_Results/01_Download/",MANIFEST$file[variable],sep=""), header=TRUE, stringsAsFactors = FALSE, sep="\t")$raw_count
}

# Copy the table for saving with UNC tags
Base -> BaseCopy

# Removes the name for writing
colnames(BaseCopy) <- gsub(".\\d\\d\\d\\d\\d\\d\\d.rsem.genes.results","", colnames(BaseCopy))

# Save the table with UNC tags
write.table(BaseCopy, file="./01_Results/BRCA.unctags.all.txt")

# Copy the table for saving with TCGA Tags
Base -> BaseCopy

#Rename Patients accordingly
colnames(BaseCopy) <- 
  SDRF$Comment..TCGA.Barcode.[
    match(colnames(BaseCopy), SDRF$Derived.Data.File)
    ]

# Creates an vector with tissue type tag
TissueTypeTag <- as.character(sapply(colnames(BaseCopy),
                          function(w) strsplit(w, "-")[[1]][4]))

# Define tissue types
table(TissueTypeTag)

# Primary Site
write.table(BaseCopy[,which(TissueTypeTag == "01A" | TissueTypeTag == "01B")],
            file="./01_Results/BRCA.PrimaryTumor.txt")

# Metastatic Site
write.table(BaseCopy[,which(TissueTypeTag == "06A")],
            file="./01_Results/BRCA.Metastatic.txt")

#Normal Tissue
write.table(BaseCopy[,which(TissueTypeTag == "11A" | TissueTypeTag == "11B")],
            file="./01_Results/BRCA.NormalTissue.txt")

# Finishes the script

