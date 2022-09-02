library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyverse)
library(REPAC)
library(edgeR)
library(limma)
# This script requires BAM files generated from simulated RNA-Seq and aligned with STAR. Simulation profiles are provided in this reposository. 
### Quantify Genes
files <- list.files("~/sims2/interaction/", pattern = "*bam$", full.names = T)

pa_quant <- featureCounts(files=files, annot.ext = "~/pa_sim_2iso_5k1k.gtf", GTF.attrType = "transcript_id",  isGTFAnnotationFile = T, strandSpecific = 0,  allowMultiOverlap = F, fraction = F, countMultiMappingReads=F, isPairedEnd = T, nthreads = 8, countReadPairs = T,  requireBothEndsMapped = F)

pa.counts <- pa_quant$counts
df <- read_tsv("~/pa_sim_2iso_5k1k.gtf", col_names = c("chr", "anno", "type", "start", "end", "score", "strand", "NA1", "metadata")) %>%
    filter(type == "transcript")
df$tid <- gsub(".+(L0.+|S0.+)\"; gene_type.+", "\\1", df$metadata)
df <- df[match(rownames(pa.counts), df$tid),]
gr <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
groups <- tibble(groups=factor(c(rep("A", 5), rep("B",5), rep("C",5), rep("D",5))))
se <- SummarizedExperiment(assays = list(counts=pa.counts),
                     rowRanges = gr,
                     colData = DataFrame(groups)
                    )

se <- sort(se)
rowRanges(se)$gene_name <- gsub("^.","", rowRanges(se)$tid)
idx <- c(names(se[strand(se) == "+"]), rev(names(se[strand(se) == "-"])))
se <- se[idx,]
rownames(se) <- unlist(tapply(rowData(se)$gene_name, rowData(se)$gene_name,
                              function(x){paste(x, sprintf("%02d", 1:length(x)), sep="_")})[unique(rowData(se)$gene_name)])
se <- se[sort(names(se))]

### Filter by expression
keep <- filterByExpr(se, group = se$groups, min.count = 10)
table(keep)
se <- se[keep,]

# Test DPU
future::plan(multisession, workers = 10)
tb.g.corr <- REPAC::fit_repac(se, group="groups", gene_name = rowData(se)$gene_name, covariates = NULL)
future::plan(sequential)


groups <- se$groups
dMat <- model.matrix(~ 0+groups)
colnames(dMat) <-  gsub("groups", "",colnames(dMat))
cMat <- makeContrasts(
    levels=colnames(dMat),
    AvsB= B - A,
    CvsD = D - C,
    g1vsg2 = ((B-A) - (D-C))
)

tb.g.corr <- fit_repac(se, design=dMat, contrasts = cMat)
write_csv(tb.g.corr$g1vsg2, "interaction_sim.csv")

int <- tb.g.corr$g1vsg2