# Statistical analysis of differential transcript usage following Salmon quantification
# https://bioconductor.riken.jp/packages/3.10/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html



#############################################################################
##   IMPORTANT: NAMES OF SAMPLES SHOULD NOT CONTAIN SPECIAL CHARS          ##
##             "_" MIGHT BE OK, "-" MESSES UP PLOTS                        ##
#############################################################################


library("tximport")
library("DRIMSeq")
library("GenomicFeatures")
library("stageR")
library("ggplot2")
library("argparse")

parser <- ArgumentParser(description = "")

parser$add_argument("--input1", type= "character")
parser$add_argument("--out_path", type= "character")
parser$add_argument("--input2", type= "character")
parser$add_argument("--alpha", type= "double")
parser$add_argument("--mge", type= "double")
parser$add_argument("--mfe", type= "double")
parser$add_argument("--mfp", type= "double")

xargs <- parser$parse_args()


suppressWarnings(dir.create(paste("results", xargs$out_path, "DRIMSeq/plots", sep="/")))

# read sample information from csv file
samps <- read.table(xargs$input1, sep=",", header=TRUE)
head(samps)

# convert the "condition" column to a factor
samps$condition <- factor(samps$condition)

# create file paths using file.path
files <- paste0(file.path("data/processed/all_samples", samps$sample_id), "/", samps$sample_id, "_salmon/quant.sf")


# set file names as sample ids
names(files) <- samps$sample_id

# import transcript quantification data with tximport
txi <- tximport(files, type="salmon", txOut=TRUE,
    countsFromAbundance="scaledTPM")

# extract transcript counts
cts <- txi$counts

# remove rows with 0 counts
cts <- cts[rowSums(cts) > 0,]

# create transcript to gene mapping
txdb <- makeTxDbFromGFF(xargs$input2, format="auto")
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")  
# count the number of transcripts per gene
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

# check if all quantified transcripts are in the database
all(rownames(cts) %in% txdf$TXNAME)

# match quantified transcripts to gene ids
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]

# check if matching is successful
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME) 

# create a data frame with gene ids, transcript ids and counts
counts <- data.frame(gene_id=txdf$GENEID,
    feature_id=txdf$TXNAME, cts, check.names = FALSE)

# create a DTU object from the counts and sample information
d <- dmDSdata(counts=counts, samples=samps)

# filter out low-expression features based on minimum number of samples and expression levels
n <- nrow(samps)
n.small <- min(aggregate(sample_id ~ condition, samps,
    function(x) length(unique(x)))$sample_id)
d <- dmFilter(d,
                min_samps_feature_expr=n.small,
                min_feature_expr=xargs$mfe,
                min_samps_feature_prop=n.small,
                min_feature_prop=xargs$mfp,
                min_samps_gene_expr=n,
                min_gene_expr=xargs$mge)

# bin gene ids based on the number of transcripts
table(table(counts(d)$gene_id))

# create the design matrix for the full model
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))

#for testing purposes for less computations
# d <- d[1:120,]

# Calculate the precision for d using the full design
set.seed(1)
system.time({
  d <- dmPrecision(d, design=design_full)
  d <- dmFit(d, design=design_full)
  d <- dmTest(d, coef="conditionem")
})


# Extract the results of the DRIMSeq analysis into res
res <- DRIMSeq::results(d)

# Extract the results of the DRIMSeq analysis at the transcript level into res.txp
res.txp <- DRIMSeq::results(d, level="feature")

# Define the no.na function to replace NAs with 1
no.na <- function(x) ifelse(is.na(x), 1, x)

# Replace NAs in res$pvalue with 1
res$pvalue <- no.na(res$pvalue)

# Replace NAs in res.txp$pvalue with 1
res.txp$pvalue <- no.na(res.txp$pvalue)

# Write res as a table to the file specified by xargs$output_TSV1 with "drimseq_genes_DTU.tsv" as suffix
write.table(res,
    file = paste("results", xargs$out_path, "DRIMSeq/genes_DTU.tsv",sep="/"),
    row.names=FALSE, na="",col.names=TRUE, sep="\t")

# Write res.txp as a table to the file specified by xargs$output_TSV2 with "drimseq_transcripts_proportions.tsv" as suffix
write.table(res.txp,
    file = paste("results", xargs$out_path, "DRIMSeq/transcripts_proportions.tsv",sep="/"),
    row.names=FALSE, na="",col.names=TRUE, sep="\t")

save(d, file = paste("results", xargs$out_path, "DRIMSeq/d.RData", sep="/"))

# Find the indices of the significant results with adjusted p-value less than xargs$alpha
signif_idx <- which(res$adj_pvalue < xargs$alpha)

# Plot the proportion for each significant result and save the plot as pdf
for (i in 1:length(signif_idx)){
    idx = signif_idx[i]
    pdf(paste("results", xargs$out_path, "DRIMSeq/plots",paste(res$gene_id[idx],"pdf",sep="."),sep="/"))
    print(plotProportions(d, res$gene_id[which(res$adj_pvalue < xargs$alpha)[i]], group_variable="condition", plot_type = "boxplot1"))
    dev.off()
}

# Prepare the data for StageR analysis
pScreen <- res$pvalue

# Define a function to shorten the gene/transcript names
strp <- function(x) substr(x,1,15)

# Shorten the gene names for pScreen
names(pScreen) <- strp(res$gene_id)

# Convert res.txp$pvalue into a matrix pConfirmation
pConfirmation <- matrix(res.txp$pvalue, ncol=1)

# Set row names of "pConfirmation" to cleaned "feature_id" from "res.txp"
rownames(pConfirmation) <- strp(res.txp$feature_id)

# Create data frame "tx2gene" with only "feature_id" and "gene_id" from "res.txp"
tx2gene <- res.txp[,c("feature_id", "gene_id")]

# Clean "feature_id" and "gene_id" in "tx2gene"
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

# Create "stageRObj" using "stageRTx" on "pScreen", "pConfirmation", and "tx2gene"
# Adjust "stageRObj" with "stageWiseAdjustment" and "dtu" method, using "stageR_alpha" from "xargs"
# Get adjusted p-values from "stageRObj" and store in "drim.padj", suppressing warnings
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                         pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=xargs$alpha)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                      onlySignificantGenes=FALSE)
})

# Write "drim.padj" to "StageR_DRIMSeq.tsv"
write.table(drim.padj,
    file = paste("results", xargs$out_path, "DRIMSeq/StageR.tsv",sep="/"),
    row.names=FALSE, na="",col.names=TRUE, sep="\t")
