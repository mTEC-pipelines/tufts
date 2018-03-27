###############################################################################
# Author: Aparna Rajpurkar
# Date: 12.15.17
#
# Required custom scripts: 
#       pipelineHead.R
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

source("funcs.R", print.eval=TRUE)
source("pipelineHead.R")


############################
### Genesets of Interest ###
############################

goi_corrplot<-c(
                "Tas2r108", "Tas2r138", "Tas2r137", "Tas2r118", 
                "Tas2r113", "Tas2r105", "Tas2r102", "Tas2r115", 
                "Tas2r107", "Tas2r104", "Tas2r117", "Tas2r123", 
                "Tas2r119", "Tas2r103", "Tas2r116", "Tas2r109"
                )

goi_heatmap<-c(
               "Tas2r108", "Tas2r138", "Tas2r137", "Tas2r118", 
               "Tas2r113", "Tas2r105", "Tas2r102", "Tas2r115", 
               "Tas2r107", "Tas2r104", "Tas2r117", "Tas2r123", 
               "Tas2r119", "Tas2r103", "Tas2r116", "Tas2r109", 
               "Trpm5", "Il17rb", "Gnat3", "Il25", "Dclk1"
               )


####################
### Make Figures ###
####################

# post normalization regressing out mito
pipeline_figures(
                 norm.dat = norm.exprs,
                 regressOut = "total_counts_MITO", 
                 sampleDat = sampleData, 
                 colorCol = c("total_counts", 
                              "total_features", 
                              "total_counts_MITO", 
                              "pct_counts_MITO"
                              ), 
                 colorGenes = c("Gnat3"), 
                 colorDiscrete = c("Genotype"),
                 titleBase = "Figures_All_Data_RegressMito_DE", 
                 goi_heatmap = goi_heatmap, 
                 goi_corrplot = goi_corrplot,
                 plotHeatGen = TRUE
                 )

# wt only
wt.sampleData = subset(sampleData, Genotype=="WT")
wt_cells = rownames(wt.sampleData)
wt.norm.exprs = norm.exprs[,wt_cells,drop=FALSE]

pipeline_figures(
                 norm.dat = wt.norm.exprs,
                 regressOut = "total_counts_MITO", 
                 sampleDat = wt.sampleData, 
                 colorCol = c("total_counts", 
                              "total_features", 
                              "total_counts_MITO", 
                              "pct_counts_MITO"
                              ), 
                 colorGenes = c("Gnat3"),
                 titleBase = "Figures_WTOnly_RegressMito_DE", 
                 goi_heatmap = goi_heatmap, 
                 goi_corrplot = goi_corrplot
                 )

# aireko only
aireko.sampleData = subset(sampleData, Genotype=="AIREKO")
aireko_cells = rownames(aireko.sampleData)
aireko.norm.exprs = norm.exprs[,aireko_cells,drop=FALSE]

pipeline_figures(
                 norm.dat = aireko.norm.exprs,
                 regressOut = "total_counts_MITO", 
                 sampleDat = aireko.sampleData, 
                 colorCol = c("total_counts", 
                              "total_features", 
                              "total_counts_MITO", 
                              "pct_counts_MITO"
                              ), 
                 colorGenes = c("Gnat3"),
                 colorDiscrete = c("Dataset"),
                 titleBase = "Figures_AIREKOOnly_RegressMito_DE", 
                 goi_heatmap = goi_heatmap, 
                 goi_corrplot = goi_corrplot
                 )
