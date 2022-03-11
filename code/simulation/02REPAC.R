library(Rsubread)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyverse)
library(REPAC)
library(edgeR)

# This script requires BAM files generated from simulated RNA-Seq and aligned with STAR. Simulation profiles are provided in this reposository. 
### Quantify Genes
files <- list.files("~/simulation/sim2iso/", pattern = "*bam$", full.names = T)

pa_quant <- featureCounts(files=files, annot.ext = "~/pa_sim_2iso_5k1k.gtf", GTF.attrType = "transcript_id",  isGTFAnnotationFile = T, strandSpecific = 1,  allowMultiOverlap = F, fraction = F, countMultiMappingReads=F, isPairedEnd = T, nthreads = 8, countReadPairs = T,  requireBothEndsMapped = F)

pa.counts <- pa_quant$counts
df <- read_table("~/pa_sim_2iso_5k1k.gtf", col_names = c("chr", "anno", "type", "start", "end", "score", "strand", "NA1", "metadata")) %>%
    filter(type == "transcript")
gr <- makeGRangesFromDataFrame(df)
gr$gene_name <- gsub("^.", "", rownames(pa.counts))
groups <- tibble(groups=factor(c(rep("A", 5), rep("B",5))))
se <- SummarizedExperiment(assays = list(counts=pa.counts),
                     rowRanges = gr,
                     colData = DataFrame(groups)
                    )

### Filter by expression
keep <- filterByExpr(se, group = se$groups, min.count = 10)
table(keep)
se <- se[keep,]

# Test DPU
future::plan(multisession, workers = 10)
tb.g.corr <- REPAC::fit_repac(se, group="groups", gene_name = rowData(se)$gene_name, covariates = NULL)
future::plan(sequential)

write_csv(tb.g.corr, file="simStats2iso.csv")
