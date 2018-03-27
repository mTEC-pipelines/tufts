###############################################################################
# Author: Aparna Rajpurkar
# Date: 12.15.17
#
# Required custom scripts: 
#       pipelineHead.R
#       MainPipelineDriver.R
# Required data files:
#       all_fc_counts.counts (count file)
#       variance_stabilized_normalization.RData (normalized data)
#       mouse_ensembl_genesymbol_conversion.txt (ENSEMBL gene symbol mapping)
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

library(limma)
library(dplyr)
library(gplots)
library(corrplot)
library(scran)
library(Rtsne)
library(RColorBrewer)
library(colorRamps)
library(reshape2)
library(pheatmap)
library(ggrepel)

# conversion file
ENSEMBL_MOUSE_CONVERSION_FILE<-"tuft_analysis_datafiles/mouse_ensembl_genesymbol_conversion.txt"
conversion<-read.delim(ENSEMBL_MOUSE_CONVERSION_FILE, header=T)

convertNames<-function(gene_list, 
                       convert_from="mgi_symbol", 
                       convert_to="ensembl_gene_id") {
    result<-conversion[match(gene_list, conversion[,convert_from]), convert_to]
    result
}

getResid<-function(norm.dat, regressOut, sampleDat) {
    regressFrmla<-"~"
    if ( length(regressOut) > 1 ) {
        concatVars<-paste(regressOut, collapse="+")
        regressFrmla<-paste0(regressFrmla, concatVars)
    } else {
        regressFrmla<-paste0(regressFrmla,regressOut)
    }
    
    regressFrmla<-as.formula(regressFrmla)

    mod<-model.matrix(regressFrmla,sampleDat)
    fit<-eBayes(lmFit(norm.dat, mod))
    resid<-residuals(fit, norm.dat)
    resid
}

getPCA<-function(exprs, scale, perform_t) {
    pca<-NULL
    if (perform_t) {
        pca<-prcomp(t(exprs), scale=scale)
    } else {
        pca<-prcomp(exprs, scale=scale)
    }
    pca
}

histCount<-function(exprs, goi, title, colorBars) {
    sum_of_goi<-as.data.frame(apply(exprs[rownames(goi),] > 0,2,sum))
    colnames(sum_of_goi)<-"total_goi_present"
    plotdatabar<-data.frame(sum_of_goi, cells=rownames(sum_of_goi))
    o<-order(plotdatabar$total_goi_present, decreasing=TRUE)
    plotdatabar$cells<-factor(plotdatabar$cells,levels=plotdatabar$cells[o])

    plot_bars<-ggplot(plotdatabar, aes(total_goi_present)) +
        geom_histogram(fill=colorBars, binwidth=1) +
        theme_classic() +
        xlab("Number of Taste Receptors Present")+
        ylab("Number of Cells")+
        ggtitle(paste("Histogram of taste receptor genes present:",title))

    plot_bars
}

ecdfplot<-function(exprs, goi, title){
    # FIXME: CLEAN THIS FUNCTION UP
    # Note: This function is rough and needs code cleaning
    # but works to produce all ecdf plots
    # it will take 10-15 mins

    ensembl_goi<-rownames(goi)

    bins_perc <- seq(0,1,0.2)
    gene_means<-as.data.frame(rowMeans(exprs))
    percentiles_all<-ecdf(gene_means[,1])
    allgenes_perc_ranks<-as.data.frame(percentiles_all(gene_means[,1]))
    rownames(allgenes_perc_ranks) <- rownames(exprs)

    genes_in_bin<-rownames(allgenes_perc_ranks)
    goi_in_bin<-genes_in_bin[genes_in_bin %in% ensembl_goi]
    cor_all_in_bin <- cor(t(exprs[genes_in_bin,]),method="spearman")

    upper_all<-upper.tri(cor_all_in_bin)
    cor_list_all<-cor_all_in_bin[upper_all]
    percentiles_cor<-ecdf(cor_list_all)
    all_cor_perc<-percentiles_cor(cor_list_all)
    cor_goi_in_bin_only<-cor_all_in_bin[goi_in_bin,goi_in_bin]
    upper_goi_in_bin<-upper.tri(cor_goi_in_bin_only)

    set_tri_na<-cor_goi_in_bin_only
    set_tri_na[upper_goi_in_bin]<-NA
    pairs_corr<-as.data.frame(as.table(set_tri_na))
    pairs_corr<-na.omit(pairs_corr)

    cor_list_goi_in_bin<-cor_goi_in_bin_only[upper_goi_in_bin]

    goi_in_bin_rank_perc<-percentiles_cor(cor_list_goi_in_bin)

    goi_in_bin.df<-cbind.data.frame(cor_list_goi_in_bin,goi_in_bin_rank_perc)
    colnames(goi_in_bin.df)<-c("Corr", "Perc")

    goi_sig<-goi_in_bin.df[(1 - goi_in_bin.df$Perc) < 0.05,,drop=FALSE]

    pairs_sig<-pairs_corr[pairs_corr$Freq %in% goi_sig$Corr,,drop=FALSE]

    pairs_sig$Var1<-as.character(convertNames(pairs_sig$Var1, 
                                              convert_from="ensembl_gene_id", 
                                              convert_to="mgi_symbol"))
    pairs_sig$Var2<-as.character(convertNames(pairs_sig$Var2, 
                                              convert_from="ensembl_gene_id", 
                                              convert_to="mgi_symbol"))

    matrix_dims <- dim(cor_goi_in_bin_only)
    pval_matrix<-as.vector(cor_goi_in_bin_only)
    pval_matrix<-1-percentiles_cor(pval_matrix)
    dim(pval_matrix)<-matrix_dims

    pval_matrix<-as.data.frame(pval_matrix)

    rownames(pval_matrix)<-as.character(convertNames(rownames(cor_goi_in_bin_only), 
                                                     convert_from="ensembl_gene_id", 
                                                     convert_to="mgi_symbol"))
    colnames(pval_matrix)<-as.character(convertNames(colnames(cor_goi_in_bin_only), 
                                                     convert_from="ensembl_gene_id", 
                                                     convert_to="mgi_symbol"))

    write.table(pval_matrix, paste0("FINAL_pval_matrix",title,".txt"),quote=FALSE)

    pdf(paste0("FINAL_",title,"_combined_ecdf_plots.pdf"))

    plot(percentiles_cor, main=paste("ECDF of all pairwise correlations,", 
                                     "Taste Receptors marked"), 
         xlim=c(-1,1), xlab="Spearman Correlation Values")

    points(goi_in_bin.df$Corr, goi_in_bin.df$Perc, pch=19, cex=1.5, col="black")
    points(goi_sig$Corr, goi_sig$Perc, pch=19, cex=1.5, col="red")
    dev.off()
}

loadingsPlot<-function(pca, genesOfInterest, title) {
    loadings<-cbind.data.frame(pca$rotation[,1:2], Genes=rownames(pca$rotation), 
                               stringsAsFactors=FALSE)

    loadings<-loadings[order(-abs(loadings$PC1),-abs(loadings$PC2)),]
    loadings<-loadings[1:50,]
    loadings$Genes<-convertNames(loadings$Genes, 
                        convert_from="ensembl_gene_id", 
                        convert_to="mgi_symbol")

    goi<-cbind.data.frame(rep("FALSE", nrow(loadings), stringsAsFactors=FALSE), 
                          stringsAsFactors=FALSE)
    
    rownames(goi)<-loadings$Genes
    goi[rownames(goi) %in% genesOfInterest[,1],]<-"TRUE"
    colnames(goi)<-"GenesOfInterest"
    loadings<-cbind.data.frame(loadings,goi)

    write.table(loadings$Genes, quote=FALSE, row.names=FALSE, col.names=FALSE, 
                paste0(title, "_loadings.txt"))

    plt<-ggplot(loadings, aes(PC1, PC2, label=Genes)) +
        geom_point(aes(color=GenesOfInterest)) +
        theme_classic() +
        geom_text_repel(aes(label=Genes,color=GenesOfInterest)) +
        ggtitle(paste(title, "loadings plot")) +
        scale_color_manual(values=c("black","red"))

    print(plt)
}

getPCA_plot_data<-function(pca, sampleDat, columnsOfSample, discreteCols, 
                           norm.genes) {
    subCols<-sampleDat[,columnsOfSample]

    if (discreteCols != "") {
        discCols<-lapply(sampleDat[,discreteCols,drop=FALSE],factor)
        plotdat<-cbind.data.frame(pca$x, subCols, discCols, t(norm.genes))
    }
    else {
        plotdat<-cbind.data.frame(pca$x, subCols, t(norm.genes))
    }
    plotdat
}

PCA_plot<-function(pcaPlotDat, xAxis, yAxis, colorCol, discreteGroup, title) {
    plt<-NULL
    if (discreteGroup == "") {
        plt<-ggplot(pcaPlotDat, aes_string(xAxis, yAxis, col=colorCol)) +
            geom_point(size=3) +
            theme_classic() +
            ggtitle(title) +
            scale_color_gradient(low = "#F0E442", high ="#009E73")
    }
    else {
        plt1<-ggplot(pcaPlotDat, aes_string(xAxis, yAxis, col=discreteGroup)) +
            geom_point(size=3, aes_string(shape=discreteGroup)) +
            theme_classic() +
            scale_color_manual(values = c("#0072B2", "#D55E00"))

        plt<-ggplot(pcaPlotDat, aes_string(xAxis, yAxis, col=colorCol)) +
            geom_point(size=3, aes_string(shape=discreteGroup)) +
            theme_classic() +
            ggtitle("Gnat3") +
            scale_color_gradient(low = "#F0E442", high ="#009E73")

        print(plt1)
    }
    print(plt)
}

cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value   
            }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
}

makeCorrPlot <- function(norm.genes) {
    correlation <- cor(t(norm.genes), method="spearman")
    diag(correlation) <- NA
	
    corrplot(correlation, method = "circle", mar = c(2,0,1,1), tl.col = "black", 
             type = "upper", col = brewer.pal(n=10, name = "RdBu"), 
             na.label = "1")
}

plotCorr<-function(norm.dat, genes) {
    colnames(genes)<-"gene_name"
    
    norm.genes<-norm.dat[rownames(genes),,drop=FALSE]
    rownames(norm.genes)<-genes$gene_name

    makeCorrPlot(norm.genes)
}

plotHeatmapGenotypes<-function(norm.dat, genes, sampleDat) {
    colnames(genes)<-"gene_name"

    norm.genes<-norm.dat[rownames(genes),,drop=FALSE]
    rownames(norm.genes)<-genes$gene_name

    heat.vals<-norm.genes - rowMeans(norm.genes)
    rownames(heat.vals)<-rownames(norm.genes)

    breaks_mid<-seq(-1,1,by=.01)
    breaks_low<-seq(min(heat.vals), -1, length.out=100)
    breaks_high<-seq(1, max(heat.vals), length.out=100)
    breaks<-unique(c(breaks_low, breaks_mid, breaks_high))
    bluewhiteyellow<-colorRampPalette(c("darkblue", "lightgrey", 
                                        "goldenrod"))(n=(length(breaks)-1))

    pheatmap(as.matrix(heat.vals), 
             color=bluewhiteyellow,
             border_color=NA,
             show_colnames=FALSE,
             cluster_rows=FALSE,
             breaks=breaks)

    pheatmap(as.matrix(heat.vals), 
            color=bluewhiteyellow,
            border_color=NA,
            show_colnames=FALSE,
            breaks=breaks)

    anno <- as.data.frame(factor(sampleDat$Genotype))
    rownames(anno) <- rownames(sampleDat)
    colnames(anno) <- "Genotype"

    # wt only clusters
    wt.cells<-rownames(subset(sampleDat, Genotype=="WT"))
    wt.heat.vals<-heat.vals[,wt.cells]
    clustercells.wt<-hclust(dist(t(as.matrix(heat.vals[,wt.cells]))));
    cluster.o.wt<-clustercells.wt$order
    wt_names_ordered<-colnames(wt.heat.vals[,cluster.o.wt])

    aireko.cells<-rownames(subset(sampleDat, Genotype=="AIREKO"))
    aireko.heat.vals<-heat.vals[,aireko.cells]
    clustercells.aireko<-hclust(dist(t(as.matrix(heat.vals[,aireko.cells]))));
    cluster.o.aireko<-clustercells.aireko$order
    aireko_names_ordered<-colnames(aireko.heat.vals[,cluster.o.aireko])

    names_ordered<-c(wt_names_ordered, aireko_names_ordered)

    heat.vals<-as.matrix(heat.vals)
    heat.vals<-heat.vals[,names_ordered]
    anno<-anno[names_ordered,,drop=FALSE]

    pheatmap(heat.vals, 
            color=bluewhiteyellow,
            border_color=NA,
            show_colnames=FALSE,
            gaps_col = length(wt.cells),
            annotation_col = anno,
            annotation_names_col = FALSE,
            cluster_cols=FALSE,
            breaks=breaks)

}

pipeline_figures<-function(norm.dat, regressOut, sampleDat, colorCol, colorGenes, 
                           colorDiscrete="", titleBase, goi_heatmap, 
                           goi_corrplot, plotHeatGen=FALSE) {
    # convert gene names to ENSEMBL
    geneNames<-rownames(norm.dat)
    ENSEMBL<-geneNames
    rownames(norm.dat)<-ENSEMBL

    # make data frames for genes of interest with both ensembl and gene symbols
    ## heatmap genes
    goi_heatmap<-data.frame(gene_name=goi_heatmap, stringsAsFactors=FALSE)
    rownames(goi_heatmap)<-convertNames(goi_heatmap$gene_name)
    ## corrplot genes
    goi_corrplot<-data.frame(gene_name=goi_corrplot, stringsasfactors=FALSE)
    rownames(goi_corrplot)<-convertNames(goi_corrplot$gene_name)
    ## genes to color onto PCA
    colorGenes<-data.frame(gene_name=colorGenes, stringsAsFactors=FALSE)
    rownames(colorGenes)<-convertNames(colorGenes$gene_name)

    # get residuals from regressing out confounders
    resid<-NULL
    if (identical(regressOut, "")) {
        resid<-norm.dat
    } else {
        resid<-getResid(norm.dat, regressOut, sampleDat)
    }

    # filter genes of interest by those that are present in data post regression
    goi_heatmap<-goi_heatmap[intersect(rownames(goi_heatmap),rownames(resid))
                             ,,drop=FALSE]
    goi_corrplot<-goi_corrplot[intersect(rownames(goi_corrplot),rownames(resid))
                           ,,drop=FALSE]
    colorGenes<-colorGenes[intersect(rownames(colorGenes), rownames(resid))
                           ,,drop=FALSE]

    # get PCA of all genes
    pca<-getPCA(resid, scale=FALSE, perform_t=TRUE)

    # get expression values for genes to color on to PCA
    colorGenesExprsValues<-resid[rownames(colorGenes),,drop=FALSE]
    rownames(colorGenesExprsValues)<-convertNames(rownames(colorGenesExprsValues),
                                                  convert_from="ensembl_gene_id",
                                                  convert_to="mgi_symbol")

    # get data to plot PCA, including color-by columns and values
    pcaplotdat<-getPCA_plot_data(pca, sampleDat, colorCol, colorDiscrete, 
                                 colorGenesExprsValues)

    ##########
    ## PDFs ##
    ##########

    pdf(paste0("FINAL_PCAs_",titleBase,".pdf"))
    # print PCAs with variables of interest mapped to them
    for (color in c(colorCol,colorGenes$gene_name)) {
        for (group in colorDiscrete) {
            PCA_plot(pcaplotdat, "PC1", "PC2", color, group, 
                     paste(titleBase, color))
        }
    }

    dev.off()

    pdf(paste0("FINAL_CORR_HEAT_HISTS",titleBase,".pdf"))

    # plot correlation values for taste receptors
    plotCorr(resid, goi_corrplot)

    # plot heatmap
    if (plotHeatGen) {
        plotHeatmapGenotypes(resid, goi_heatmap, sampleDat)
    }

    # plot hustograms of taste receptor genes
    print(histCount(resid, goi_heatmap, titleBase, "#D55E00"))
    print(histCount(resid, goi_heatmap, titleBase, "#0072B2"))

    dev.off()
    
    # plot ECDFs
    ecdfplot(resid, goi_corrplot, titleBase)
}
