library(SingleCellExperiment)
library(scater)

MIN_READS = 100000
MIN_FEATURES = 750
MAX_MITO_PERC = 10

#############
# arguments #
#############

directory <- "tuft_analysis_datafiles/"

aireko_file <- "aireko_counts.txt"

wt_file <- "all_wt_counts.txt"
pdf_name <- "all_counts_histograms.pdf"

wt.expression.data <- read.delim(paste(directory, wt_file, sep = ""), header=T)
aireko.expression.data <- read.delim(paste(directory, aireko_file, sep = ""), header=T)

# create batches for wt
numCells.wt<-ncol(wt.expression.data)
batch <-c(rep("1",numCells.wt))
Feb_ids<-grep("tuft", colnames(wt.expression.data))
batch[Feb_ids]<-"2"
genotype <- c(rep("WT",numCells.wt))
cell_info_wt <- data.frame(batch=batch, genotype=genotype) # Must be data frame object
rownames(cell_info_wt) <- colnames(wt.expression.data)

# create batches for aire KO
numCells.aireko <- ncol(aireko.expression.data)
batch <- c(rep("2", numCells.aireko))
genotype <- c(rep("AireKO",numCells.aireko))
cell_info_aireko <- data.frame(batch=batch, genotype=genotype)
rownames(cell_info_aireko) <- colnames(aireko.expression.data)

# Make one file containing all count data
cell_info <- rbind(cell_info_wt, cell_info_aireko)

expression.data <- merge(wt.expression.data, aireko.expression.data, by = "row.names", all=TRUE)

rownames(expression.data) <- expression.data$Row.names

expression.data$Row.names <- NULL

##########################################################################

# deal with the ENSEMBL gene name issue
ensembl_mouse<-rownames(expression.data)[grepl("ENSMUSG", rownames(expression.data))]   
rownames(expression.data)<-sub("(ENSMUSG\\d+|ERCC.+)\\.*\\d*", "\\1", rownames(expression.data), perl=TRUE)        

expression.data_int <- as.matrix(expression.data)

# create SCESet to hold all data
mtecdev<-SingleCellExperiment(assays= list(counts=expression.data_int), colData=cell_info)

keep_feature <- rowSums(counts(mtecdev) > 0) > 0
mtecdev <- mtecdev[keep_feature, ]

# set ERCC and mito
isSpike(mtecdev, "ERCC") <- grepl("ERCC", rownames(mtecdev))
isSpike(mtecdev, "MT") <- rownames(mtecdev) %in%
  c('ENSMUSG00000064336', 'ENSMUSG00000064337', 'ENSMUSG00000064338',
  'ENSMUSG00000064339', 'ENSMUSG00000064340', 'ENSMUSG00000064341',
  'ENSMUSG00000064342', 'ENSMUSG00000064343', 'ENSMUSG00000064344',
  'ENSMUSG00000064345', 'ENSMUSG00000064346', 'ENSMUSG00000064347',
  'ENSMUSG00000064348', 'ENSMUSG00000064349', 'ENSMUSG00000064350',
  'ENSMUSG00000064351', 'ENSMUSG00000064352', 'ENSMUSG00000064353',
  'ENSMUSG00000064354', 'ENSMUSG00000064355', 'ENSMUSG00000064356',
  'ENSMUSG00000064357', 'ENSMUSG00000064358', 'ENSMUSG00000064359',
  'ENSMUSG00000064360', 'ENSMUSG00000064361', 'ENSMUSG00000065947',
  'ENSMUSG00000064363', 'ENSMUSG00000064364', 'ENSMUSG00000064365',
  'ENSMUSG00000064366', 'ENSMUSG00000064367', 'ENSMUSG00000064368',
  'ENSMUSG00000064369', 'ENSMUSG00000064370', 'ENSMUSG00000064371',
  'ENSMUSG00000064372')


# calculate basic QC, setting controls as ERCCs 
mtecdev <- calculateQCMetrics (
  mtecdev,
  feature_controls = list(
    ERCC = isSpike(mtecdev, "ERCC"),
    MT = isSpike(mtecdev, "MT")
  )
)

pdf(paste(directory, pdf_name, sep = ""))



# Make plotting df for ggplot
df_to_plot <- data.frame(total_features=mtecdev$total_features,
  total_counts = mtecdev$total_counts, mito_pct= mtecdev$pct_counts_MT,
  ercc_pct = mtecdev$pct_counts_ERCC, Batch = mtecdev$batch,
  Genotype = mtecdev$genotype)

# Histogram of total counts
counts_plot <- ggplot(df_to_plot, aes(total_counts))
counts_plot + geom_histogram(aes(fill = Genotype), bins = 100, col = "black") + theme_classic() +
scale_color_manual(values = c("#D55E00", "#0072B2")) +
labs(x = "Total Counts") + geom_vline(xintercept = MIN_FEATURES, color = "black")

# Histogram of total genes
features_plot <- ggplot(df_to_plot, aes(total_features))
features_plot + geom_histogram(aes(fill = Genotype), bins = 100, col = "black") + theme_classic() +
scale_color_manual(values = c("#D55E00", "#0072B2")) +
labs(x = "Total Features") + geom_vline(xintercept = MIN_FEATURES, color = "black")

# Scatter plot of genes vs percent mitochondria
mito_plot <- ggplot(df_to_plot, aes(total_features, mito_pct))
mito_plot + geom_point(aes(shape=Batch, color=Genotype), size = 3) +
scale_color_manual(values = c("#D55E00", "#0072B2")) + theme_classic() +
labs(x = "Total Features", y = "Mitochondrial Read Percent") +
geom_hline(yintercept = 10, color = "black")

dev.off()

