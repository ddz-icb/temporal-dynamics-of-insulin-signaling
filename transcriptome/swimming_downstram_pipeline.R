#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#install.packages("RMariaDB")

library(tximport)
#library(DRIMSeq)
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(GenomicFeatures)

# please change the path to your own file structure
sample_path = "C:\\your_path_to_working_directory\\"

# experiment setup
samps <- read.csv(paste(c(sample_path, "samples2.csv"), sep = "",collapse=""),sep=",",header = TRUE)
samps$condition <- factor(samps$condition)
table(samps$condition)

# load abundance files
files <- file.path(sample_path,"samples", samps$sample_id, "abundance.tsv", fsep="\\")
names(files) <- samps$sample_id
#head(files)   
txi <- tximport(files, type="kallisto", txOut=TRUE,
                countsFromAbundance="scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]
df <- as.data.frame(cts)
row.names(df) <-sub('\\.[0-9]*$', '',rownames(cts)) #gtf to txdf is missing transcript version id
cts2 <- as.matrix(df)

# annotation
#txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
gtf <- "Homo_sapiens.GRCh38.110.gtf"
txdb3 <- makeTxDbFromGFF(gtf,
                format="gtf",
                dataSource="Ensembl",
                organism="Homo sapiens",
                taxonomyId=NA,
                circ_seqs=NULL)
txdf <- select(txdb3, keys(txdb3, "GENEID"), "TXNAME", "GENEID")
# txdb2 <- makeTxDbFromEnsembl(organism="Homo sapiens",
#                              release=108,
#                              circ_seqs=NULL,
#                              server="ensembldb.ensembl.org",
#                              username="anonymous", password=NULL, port=0L,
#                              tx_attrib=NULL)
#txdf <- select(txdb2, keys(txdb2, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID) # count number of transcripts
txdf$ntx <- tab[match(txdf$GENEID, names(tab))] # add number of transcripts
#txdf2 <- txdf[match(rownames(cts2),txdf$TXNAME),]
txdf2 <- as.data.frame(txdf)

#DRIMSeq
counts <- merge(txdf2, cts2, by.x ="TXNAME", by.y = 0)
library(dplyr)
counts <- counts %>% rename(feature_id = TXNAME,gene_id =  GENEID)
library(DRIMSeq)
d <- dmDSdata(counts=counts, samples=samps)

#Filtering
n <- 12
n.small <- 6
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)

design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)
set.seed(1)
system.time({
  d <- dmPrecision(d, design=design_full)
  d <- dmFit(d, design=design_full)
  d <- dmTest(d, coef="conditionInsulin")
})

# statitic
res <- DRIMSeq::results(d) # test any differential transcript usage within the gene
res.txp <- DRIMSeq::results(d, level="feature")
# test single p-value per transcript, which tests
# whether the proportions for this transcript changed within the gene

# replace na p-vlues with 1
no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)

# plot one gene
idx <- which(res$adj_pvalue < 0.05)[1]
res[idx,]
#FASN = ENSG00000169710
# Smarce1 ENSG00000073584
# gmeb1 ENSG00000162419
# rgl4 ENSG00000159496
idx <- which(res$gene_id == "ENSG00000073584")
res[idx,]

#plotProportions(d, res$gene_id[idx], "condition",plot_type="ribbonplot")
#trace(plotProportions, edit=TRUE)
#png(file="C:\\Users\\smajda\\Documents\\PROJEKT_phospho_1\\results\\figures\\saving_plot_Smarce1.png",
#    width=800, height=850)
dm_plotProportions(d, res$gene_id[idx], "condition",plot_type="boxplot1")
#dev.off

# final annotation
anno <- read.csv("biomart_gene_2_go_name_short.txt",sep="\t",header = TRUE)
anno <- as.data.frame(anno)
merged_df <- merge(res, anno, by.x ="gene_id", by.y = "Gene.stable.ID", all.x = TRUE)
write.table(merged_df,"all_DTU_genes_anno.tsv", row.names = FALSE, sep = "\t")

# for transcripts
merged_df <- merge(res.txp, anno, by.x ="gene_id", by.y = "Gene.stable.ID", all.x = TRUE)
write.table(merged_df,"transcripts.tsv", row.names = FALSE, sep = "\t")

#for next R script:
anno <- read.csv("merged_table2_biotype.csv",header = TRUE)
anno <- as.data.frame(anno)
anno$ensembl_transcript_id_version <- sub('\\.[0-9]*$', '',ensembl_transcript_id_version)
anno <- anno[ , ! names(anno) %in% c("pvalue_adjusted")]
ext_table_short <- res.txp[c("feature_id","adj_pvalue")]
merged_df <- merge(ext_table_short, anno, by.x ="feature_id", by.y = "ensembl_transcript_id_version", all.x = TRUE)
write.table(merged_df,"filtered_transcripts.tsv", row.names = FALSE, sep = "\t")

# for volcano
anno <- read.csv("biomart_gene_2_go_name_short.txt",sep="\t",header = TRUE)
anno <- as.data.frame(anno)
merged_df <- merge(res, anno, by.x ="gene_id", by.y = "Gene.stable.ID", all.x = TRUE)
anno2 <- read.csv("annotation_for_volcano.csv",header = TRUE)
anno2 <- as.data.frame(anno2)
merged_df2 <- merge(merged_df, anno2, by.x ="Gene.name", by.y = "hgnc_symbol", all.x = TRUE)
write.table(merged_df2,"all_DTU_genes_anno_log2FC.tsv", row.names = FALSE, sep = "\t")

