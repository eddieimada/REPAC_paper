library(tidyverse)
library(SummarizedExperiment)
library(annotatr)
library(GenomicRanges)
library(edgeR)
library(furrr)
library(compositions)
library(car)
library(recount3)
library(REPAC)

# Select project
mouse_projects <- available_projects(organism = "mouse")
proj_info <- subset(
    mouse_projects,
    project == "SRP048707" & project_type == "data_sources"
)

# create PA RSE
df <-read_tsv("data/QAPA_up_mm10_sorted.bed", col_names = F)
df <- df[,c(1,2,3,6,7,8)]
names(df) <- c("chr", "start", "end", "strand", "gene_name", "annotation")
gr <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
gr <- gr[!duplicated(gr)]

future::plan(multisession, workers = 10)
system.time(
se <- create_pa_rse(organism = "mouse", project="SRP048707", annotation=gr)
)
future::plan(sequential)
colnames(se) <- se$external_id
se$sra.sample_title <- gsub(".$", "", se$sra.sample_title)
# Select untreated samples
se <- se[,grepl("WT", se$sra.sample_title)]
write_rds(se, "objs/SRP048707_APA_qapa.rds")

# Filter low counts
pa.counts <- assays(se)$counts
### Set groups
groups <- se$sra.sample_title
pheno <- groups
### Filter by expression
#load gene level counts from salmon
load("bcell_salmon.rda")
dge <- DGEList(mat, genes = fmat, group = se$sra.sample_title)

keep <- filterByExpr(dge, group = se$sra.sample_title, min.count= 30, min.total.count=10)
table(keep)

gnsExp <-dge$genes$external_gene_name[keep]
keep2 <- gsub("_.+", "", rownames(pa.counts)) %in% gnsExp
keep <- filterByExpr(pa.counts, group = groups, min.count = 10)
table(keep & keep2)
pa.counts <- pa.counts[keep & keep2,]
se <- se[keep & keep2,]    
se <- se[!duplicated(rowRanges(se)),]

# keep only 3'UTR sites
se <- se[rowRanges(se)$annotation == "3UTR",]

# Test for DPU
future::plan(multisession, workers = 10)
system.time(
tb.g.corr <- fit_repac(se, group="sra.sample_title", covariates = NULL)
)
future::plan(sequential)

tb.g.corr <- tb.g.corr %>%
    arrange(abs(p.val)) %>%
    filter(!duplicated(Ref))

write_csv(tb.g.corr, file="bcell.csv")



tb.ord <- tb.g.corr %>%
    filter(abs(cFC) >= 0.25 & adj.p.val <= 0.05 ) %>%
    arrange(abs(adj.p.val))
