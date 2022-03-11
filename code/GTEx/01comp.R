library(tidyverse)
library(SummarizedExperiment)
library(GenomicRanges)
library(edgeR)
library(furrr)
library(compositions)
library(car)
library(recount3)
library(megadepth)
library(REPAC)
human_samples <- available_samples(organism = "human")
gtex.samples.b <- human_samples[human_samples$project == "BRAIN",]$external_id
gtex.samples.t <- human_samples[human_samples$project == "TESTIS",]$external_id
meta.b <- locate_url("BRAIN", project_home = "data_sources/gtex", type="metadata", sample = gtex.samples.b)
meta.b <- file_retrieve(meta.b)
meta.b <- read_metadata(meta.b)
meta.b <- meta.b[meta.b$gtex.smtsd == "Brain - Cortex",]

meta.t <- locate_url("TESTIS", project_home = "data_sources/gtex", type="metadata", sample = gtex.samples.t)
meta.t <- file_retrieve(meta.t)
meta.t <- read_metadata(meta.t)

meta <- bind_rows(meta.b, meta.t)
keep <- map_lgl(meta, ~ all(!is.na(.x)))
meta <- meta[,keep]
meta$BigWigURL <- urls
isRNA <- meta$gtex.smnabtcht == "RNA Extraction from Paxgene-derived Lysate Plate Based"
isGood <- meta$gtex.smafrze == "RNASEQ"

meta <- meta[isRNA & isGood,]
set.seed(1)
samp10 <- c(
sample(meta$external_id[meta$gtex.smtsd == "Brain - Cortex"], 10),
sample(meta$external_id[meta$gtex.smtsd == "Testis"], 10)
)

samp100 <- c(
    sample(meta$external_id[meta$gtex.smtsd == "Brain - Cortex"], 50),
    sample(meta$external_id[meta$gtex.smtsd == "Testis"], 50)
)

samp200 <- c(
    sample(meta$external_id[meta$gtex.smtsd == "Brain - Cortex"], 100),
    sample(meta$external_id[meta$gtex.smtsd == "Testis"], 100)
)

samp400 <- c(
    sample(meta$external_id[meta$gtex.smtsd == "Brain - Cortex"], 200),
    sample(meta$external_id[meta$gtex.smtsd == "Testis"], 200)
)


future::plan(multisession, workers = 8)
system.time(
se10 <- create_pa_rse(organism = "human", project=c("BRAIN", "TESTIS"),
                    annotation="data/QAPA_up_hg38.bed",
                    sample_id = samp10,
                    bed_cols= c("seqnames", "start", "end", "ids", 
                                "score", "strand", "SYMBOL", "annotation",
                                "isTE", "trash"))
)

system.time(
    se100 <- create_pa_rse(organism = "human", project=c("BRAIN", "TESTIS"),
                          annotation="data/QAPA_up_hg38.bed",
                          sample_id = samp100,
                          bed_cols= c("seqnames", "start", "end", "ids", 
                                      "score", "strand", "SYMBOL", "annotation",
                                      "isTE", "trash"))
)

system.time(
    se200 <- create_pa_rse(organism = "human", project=c("BRAIN", "TESTIS"),
                           annotation="data/QAPA_up_hg38.bed",
                           sample_id = samp200,
                           bed_cols= c("seqnames", "start", "end", "ids", 
                                       "score", "strand", "SYMBOL", "annotation",
                                       "isTE", "trash"))
)

system.time(
    se400 <- create_pa_rse(organism = "human", project=c("BRAIN", "TESTIS"),
                           annotation="data/QAPA_up_hg38.bed",
                           sample_id = samp400,
                           bed_cols= c("seqnames", "start", "end", "ids", 
                                       "score", "strand", "SYMBOL", "annotation",
                                       "isTE", "trash"))
)
future::plan(sequential)

save(se10, se100, se200, se400, file="objs/gtex_batches.rda")
load("objs/gtex_batches.rda")

walk(c("se10", "se100", "se200", "se400"), function(id) {
    se <- get(id)
colnames(se) <- se$external_id
se <- sort(se)
idx <- c(names(se[strand(se) == "+"]), rev(names(se[strand(se) == "-"])))
se <- se[idx,]
ids <- tibble(symbol=rowData(se)$SYMBOL) %>% group_by(symbol) %>%
    mutate(ids=paste(symbol, sprintf("%02d", 1:n()), sep="_"))
ids <- ids$ids
rowData(se)$ids <- ids
rownames(se) <- ids
se <- se[sort(names(se)),]

se <- se[,colSums(assays(se)$counts) > 0]
pa.counts <- assays(se)$counts
### set groups
groups <- se$gtex.smts
pheno <- groups
### Filter by expression
keep <- filterByExpr(pa.counts, min.count = 10)
table(keep)
pa.counts <- pa.counts[keep,]
se <- se[keep,]    
se <- se[!duplicated(rowRanges(se)),]
se <- se[rowRanges(se)$annotation == "3UTR",]

future::plan(multisession, workers = 8)
system.time(
tb.g.corr <- fit_corepad(se, group="gtex.smts", covariates = NULL)
)
future::plan(sequential)

write_csv(tb.g.corr, file=paste0("text/BT_", id,".csv"))

})

tb.g.corr <- read_csv("text/BT_se400.csv")
tb.g.corr <- tb.g.corr %>%
    arrange(abs(p.val)) %>%
    filter(!duplicated(Ref))

tb.ord <- tb.g.corr %>%
    filter(abs(cFC) >= 0.25 & adj.p.val <= 0.05 ) %>%
    arrange(-abs(cFC))
