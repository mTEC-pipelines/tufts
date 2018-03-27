###############################################################################
# Author: Aparna Rajpurkar
# Date: 2.15.17
#
# Required custom scripts: 
#       MainPipelineDriver.R
#       funcs.R
# Required data files:
# Data files should be contained within a subdirectory named:
#       tuft_analysis_datafiles
#
# Purpose
# 
# These scripts perform analysis for single cell analysis of Tuft cells.
# No guarantees made. These scripts were written for personal use.
# Contact: arrajpur@stanford.edu with questions.
# 
# Usage: Rscript MainPipelineDriver.R
# Modify the driver file for different plots
###############################################################################

#source("funcs.R", print.eval=TRUE)
library(scran) #needed for SingleCellExperiment
library(scater) #needed for calculateQCMetrics

ENSEMBL_MOUSE_CONVERSION_FILE<-"tuft_analysis_datafiles/mouse_ensembl_genesymbol_conversion.txt"
NORMALIZED_DATA_FILE<-"tuft_analysis_datafiles/POR_all_data.rda"
WT_COUNTS<-"tuft_analysis_datafiles/all_wt_counts.txt"
AIREKO_COUNTS<-"tuft_analysis_datafiles/aireko_counts.txt"
BATCH_INFO<-"tuft_analysis_datafiles/samples_all.txt"

conversion<-read.delim(ENSEMBL_MOUSE_CONVERSION_FILE, header=T)

convertNames<-function(gene_list, convert_from="mgi_symbol", convert_to="ensembl_gene_id") {
    result<-conversion[match(gene_list, conversion[,convert_from]), convert_to]
    result
}

# get normalized data
pos_control<-"kw_tuft_142"
vars<-load(NORMALIZED_DATA_FILE)
norm.exprs<-normalized
norm.exprs<-norm.exprs[,!(colnames(norm.exprs)==pos_control)]
norm.exprs<-norm.exprs[!(grepl("ERCC", rownames(norm.exprs))),,drop=FALSE]

# get raw counts for basic metadata
wt_file <- WT_COUNTS
aire_file <- AIREKO_COUNTS

controls<-c("kw_tuft_141", "kw_tuft_142", "kw_tuft_143", "kw_tuft_144", "June29_KW.96.featureCounts")

DATA_WT <- read.csv(wt_file, header = T, row.names = 1, sep="\t")

DATA_AIRE <- read.csv(aire_file, header = T, row.names = 1, sep="\t")
DATA <- cbind.data.frame(DATA_WT, DATA_AIRE[rownames(DATA_WT),,drop=FALSE])

raw.counts<-DATA[,!(colnames(DATA) %in% controls),drop=FALSE]
raw.counts<-raw.counts[!(grepl("ERCC", rownames(raw.counts))),,drop=FALSE]

# get total number of cells from raw counts
numCells <- ncol(raw.counts)

# set batches
batch_info<-read.delim(BATCH_INFO, header=T, row.names=1)
batch<-batch_info[colnames(raw.counts),,drop=FALSE]

##########################################################################

# deal with the ENSEMBL gene name issue
ensembl_mouse<-rownames(raw.counts)[grepl("ENSMUSG", rownames(raw.counts))]   
rownames(raw.counts)<-sub("(ENSMUSG\\d+|ERCC.+)\\.*\\d*", "\\1", rownames(raw.counts), perl=TRUE)        

geneNamesForMito<-convertNames(rownames(raw.counts), convert_from="ensembl_gene_id", convert_to="mgi_symbol")

mito_bool<-grepl("mt", geneNamesForMito)

# get ercc names from raw raw.counts
erccs_raw<-rownames(raw.counts[grepl("ERCC",rownames(raw.counts)),])
mito_raw<-rownames(raw.counts[mito_bool,])

# set phenotype data (describe cells)
pd<-data.frame(batch=batch)
rownames(pd)<-colnames(raw.counts)
colnames(pd)<-colnames(batch)

fd<-data.frame(rownames(raw.counts))
rownames(fd)<-rownames(raw.counts)
colnames(fd)<-"ensembl"

# create SingleCellExperiment to hold all data
mtecdev<-SingleCellExperiment( 
                           list(counts=as.matrix(raw.counts)),
                           colData=pd,
                           rowData=fd
                           )

# calculate basic QC, setting controls as ERCCs 
mtecdev<-calculateQCMetrics(
                            mtecdev,
                            feature_controls=list(ERCC=erccs_raw, MITO=mito_raw)
                            )


# get phenotype data / sample metadata
sampleData<-cbind.data.frame(colData(mtecdev))
rownames(sampleData)<-colnames(mtecdev)

norm.cells<-colnames(norm.exprs)

sampleData<-sampleData[norm.cells,]
