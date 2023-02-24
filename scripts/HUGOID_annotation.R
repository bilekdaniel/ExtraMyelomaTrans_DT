library("AnnotationDbi")
library("org.Hs.eg.db")
library("argparse")
library("tibble")

parser <- ArgumentParser(description = "")

parser$add_argument("--input1", type= "character")
parser$add_argument("--input2", type= "character")
parser$add_argument("--input3", type= "character")
parser$add_argument("--out_path", type= "character")

xargs <- parser$parse_args()

#function that removes .xx in ENSG ID
strp <- function(x) substr(x,1,15)

########################################
##              genes_DTU             ##
########################################

#load file
ANNOT <- read.table(xargs$input1, header = TRUE, sep="\t")
ANNOT_strip <- strp(ANNOT$gene_id)
ANNOT$gene_id_str <- ANNOT_strip

#translation from ENSG to Hugo symbol
HugoID <- mapIds(org.Hs.eg.db, keys = ANNOT$gene_id_str, keytype = "ENSEMBL", column = "SYMBOL")
ANNOT$HugoID <- HugoID

#reorder columns for more clarity
col_order <- c("gene_id", "gene_id_str", "HugoID", "lr", "df", "pvalue", "adj_pvalue")
ANNOT <- ANNOT[, col_order]

#save the file
write.csv(ANNOT, file = paste("results", xargs$out_path, "DRIMSeq/genes_DTU.csv", sep = "/"))

########################################
##              PROP                  ##
########################################

ANNOT <- read.table(xargs$input2, header = TRUE, sep="\t")
ANNOT_strip <- strp(ANNOT$gene_id)
ANNOT$gene_id_str <- ANNOT_strip

HugoID <- mapIds(org.Hs.eg.db, keys = ANNOT$gene_id_str, keytype = "ENSEMBL", column = "SYMBOL")
ANNOT$HugoID <- HugoID

col_order <- c("gene_id", "gene_id_str", "HugoID", "feature_id", "lr", "df", "pvalue", "adj_pvalue")
ANNOT <- ANNOT[, col_order]

write.csv(ANNOT, file = paste("results", xargs$out_path, "DRIMSeq/transcripts_proportions.csv", sep = "/"), row.names = FALSE)

########################################
##              StageR                ##
########################################

ANNOT <- read.table(xargs$input3, header = TRUE, sep="\t")
HugoID <- mapIds(org.Hs.eg.db, keys = ANNOT$geneID, keytype = "ENSEMBL", column = "SYMBOL")
ANNOT$HugoID <- HugoID

col_order <- c("geneID", "HugoID", "txID", "gene", "transcript")
ANNOT <- ANNOT[, col_order]

write.csv(ANNOT, file = paste("results", xargs$out_path, "DRIMSeq/StageR.csv", sep = "/"))
