library(rnaseqGene)

dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
csvfile <- file.path(dir, "sample_table.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata

coldata$names <- coldata$Run
dir <- "/home/islam/airway"
coldata$files <- file.path(dir, "quants", paste0(coldata$names,"_quant"), "quant.sf")
file.exists(coldata$files)

# makeLinkedTxome
# check the latest genomes versions
hashfile <- file.path(system.file("extdata",package="tximeta"),"hashtable.csv")
hashtable <- read.csv(hashfile,stringsAsFactors=FALSE)
hashtable[,"organism"]

se <- tximeta(coldata,useHub=FALSE)

gse <- summarizeToGene(se)

assay(gse)
assayNames(gse)
head(assay(gse), 3)

colSums(assay(gse))

rowRanges(gse)

seqinfo(rowRanges(gse))

colData(gse)

gse$dex <- factor(gse$dex, levels=c("untrt","trt"))
gse$cell <- factor(gse$cell, unique(gse$cell))

levels(gse$cell)

dds <- DESeqDataSet(gse, design = ~ cell + dex)

keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)

log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)



dds <- DESeq(dds)
res <- results(dds)
res
summary(res)
mcols(res, use.names = TRUE)


library("apeglm")
resultsNames(dds)

#results(dds, contrast = c("cell", "N061011", "N61311"))

res <- results(dds, contrast=c("dex","trt","untrt"))
summary(res)

DESeq2::plotMA(res)

res <- lfcShrink(dds, coef="dex_trt_vs_untrt", type="apeglm")
DESeq2::plotMA(res, ylim = c(-5, 5))


res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)


resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)






