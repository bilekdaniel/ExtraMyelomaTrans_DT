
# Differential expression analysis based on: 
# https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#summarizedexperiment

library("tximeta")
library("DESeq2")
library("dplyr")
library("ggplot2")
library("AnnotationDbi")
library("genefilter")
library("pheatmap")
library("ggrepel")
library("biomaRt")
library("DEGreport")
library("apeglm")


# Parse arguments from a command line
args = commandArgs(trailingOnly=TRUE)


# Read metadata table
sample_metadata <- read.table(args[3], sep=',', header=TRUE, stringsAsFactors = FALSE)


# Ad $files columns with path to the quant.sf files (result from salmon)
sample_metadata$files <- paste0(args[1],"/", sample_metadata$names, "/", sample_metadata$names, "_salmon/quant.sf")
# sample_metadata$files <- paste0("data/processed/second_run","/", sample_metadata$names, "/", sample_metadata$names, "_salmon/quant.sf")




suppressPackageStartupMessages(library(tximeta))
makeLinkedTxome(indexDir="/scratch/project/open-22-59/gencode/salmon_index",
                source="GENCODE",
                organism="Homo sapiens",
                release="13",
                genome="GRCh38",
                fasta="/scratch/project/open-22-59/gencode/gencode.v36.transcripts.fa",
                gtf="/scratch/project/open-22-59/gencode/gencode.v36.annotation.gtf",
                write=TRUE)


# Parse salmon result with tximeta
se <- tximeta(sample_metadata)

# summarize transcript counts to genes
gse <- summarizeToGene(se)
row.names(gse) = substr(rownames(gse), 1, 15)


# Construction of DESeqDataSet object
dds <- DESeqDataSet(gse, design = ~ specimen)


# Remove rows of the DESeqDataSet that have no counts, or only a single count across all samples
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]



# EXPLORATORY ANALYSIS AND VISUALIZATION
########################################

# counts normalization
rld <- rlog(dds, blind = FALSE)

# principal component analysis

pcaData <- plotPCA(rld, intgroup = c( "specimen", "patient"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(paste0("results/", args[2], "/deseq2/all/pca_rld.pdf"))
ggplot(pcaData, aes(x = PC1, y = PC2, shape = specimen, color = patient)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with rld data")
dev.off()


# Differential expression analysis
##################################


# differential expression pipeline
dds <- DESeq(dds, test="LRT", reduced=~1)



# create a result table  with False discovery rate 0.05

res <- results(dds, alpha = 0.05)
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef=2, type="apeglm")
res





counts_ = counts(dds, normalized=T)
write.csv(counts_, file=paste0("results/", args[2], "/deseq2/all/norm_counts_all.csv"))

# filter our hits with a baseMean lower than 10
res = subset(res, res$baseMean > 10)



mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name", "description"),
  filter="ensembl_gene_id",
  values=rownames(res),
  uniqueRows=TRUE)


#############################
############################
# genes with no annotation
problematic_ids = setdiff(rownames(res), annotLookup$ensembl_gene_id)

fileConn <- file( paste0("results/", args[2], "/deseq2/all/problematic_ids.out"))    
writeLines(problematic_ids, fileConn)    
close(fileConn)

res <- subset(res, !(rownames(res) %in% problematic_ids))
ens = rownames(res)

annotLookup <- data.frame(
  ens[match(annotLookup$ensembl_gene_id, ens)],
  annotLookup)

colnames(annotLookup) <- c( 
      "original_id", 
      c("ensembl_gene_id", "gene_biotype", "external_gene_name", "description")) 

res$symbol <- annotLookup$external_gene_name
res$biotype <- annotLookup$gene_biotype
res$description <- annotLookup$description


######################################
#####################################


# order according to padj
resOrdered <- res[order(res$pvalue),]
resOrderedDF <- as.data.frame(resOrdered)
# filtResOrderedDF = filter(resOrderedDF, padj < 0.05)
write.csv(resOrderedDF, file = paste0("results/", args[2], "/deseq2/all/deseq_result_all.csv"))





# plot heatmap with all significant genes
pdf(paste0("results/", args[2], "/deseq2/all/significant_genes_heatmap.pdf"),  height=14)
resSig <- subset(res, padj < 0.01)
resSig <- subset(resSig, biotype != "IG_V_gene")
resSig <- subset(resSig, biotype != "IG_C_gene")
select <- rownames(resSig[ order(resSig$padj), ], )
mat = assay(rld)[select,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[,c("specimen", "paired")])
pheatmap(mat, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=anno, show_colnames=TRUE)
dev.off()


