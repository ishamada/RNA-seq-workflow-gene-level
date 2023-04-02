# BiocManager::install("rnaseqGene")
library(rnaseqGene)
library("airway")
dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))

csvfile <- file.path(dir, "sample_table.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata

coldata <- coldata[1:2,]
coldata$names <- coldata$Run
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf.gz")
file.exists(coldata$files)

library("tximeta")
se <- tximeta(coldata)
dim(se)
head(rownames(se))

gse <- summarizeToGene(se)


rowRanges(gse)

colData(gse)

assayNames(gse)
head(assay(gse), 3)

gse$dex

gse$dex <- factor(gse$dex, levels=c("untrt","trt"))

levels(gse$dex)

gse$dex %<>% relevel("untrt")
gse$dex

gse$dex <- relevel(gse$dex, "untrt")


dds <- DESeqDataSet(gse, design = ~ dex)


