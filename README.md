# TNBC-Targeting

This GitHub Repository contains all the scripts and working examples used in the following publication:

Guanylate-binding protein-1 is a potential new therapeutic target for triple negative breast cancer
(submitted to BMC Cancer)

Authors:
Melissa Quintero, Douglas Adamoski, Larissa Menezes dos Reis, Carolline Fernanda Rodrigues Ascenção, Kaliandra de Almeida Gonçalves, Marília Meira Dias, Marcelo Falsarella Carazzolle, Sandra M. G. Dias

Details:

R-scripts (inside R-scripts folder)

Scripts should be executed in numerical order (as below), since each one create files for the next.
Seven folders will be created, containing results from that step. Don't rename them, as the next script will look for files in the folder with that name.

*Script01_GDC-Download.R:
Download TCGA files from GDC data portal

*Script02_TNBC_Assignment.R
Assign breast marker status using RNA levels and select the best mixture model comparing to IHC data

*Script03_TCGA_DESeq2-Analysis.R
Perform differential expression analysis using DESeq2 and TCGA data, comparing Triple Negative Breast Cancer (TNBC) against Non-TNBC

*Script04_TCGA_PathwayEnrichment.R
Perform pathway enrichments using DE genes from TCGA cohort

*Script05_Cell_Line_DESeq2_PathwayEnrichment.R
perform DE genes analysis and pathway enrichments using RNA-seqs from cell lines

*Script06_PubMed_Other_Graphs.R
Perform access information about genes in PubMed, canSAR database druggability and plot some graphs from experimental data

*Script07_OtherMicroArray_Databases.R
Access reprodutibiliy of TCGA breast cancer cohort in external microarray datasets

*ExternalFiles:
This folder contains all files received outside R (this includes qPCR results, cellular proliferation assays, cell line classifications, RSEM values from in-house and GEO processed RNA-seqs and canSAR database outputs)


