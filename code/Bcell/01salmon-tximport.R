library(tximport)
library(data.table)
library(Biobase)
library(biomaRt)

# This script requires salmon quantification. Please refer to the manuscript on how these files were generated.
### Load transcript to gene table
trsnc <- read_tsv("~/Downloads/bcell2/WT_CTRL1_cds/quant.sf")[,1, drop = T]
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl",version = "104")
feat <- getBM(attributes=c("ensembl_transcript_id_version",
                           "ensembl_gene_id_version",
                           "external_gene_name",
                           "description"),
              filters = "ensembl_transcript_id_version",
              values = trsnc, #keys values for query
              mart = ensembl)
feat <- feat[!duplicated(feat$ensembl_transcript_id_version),]
tx2gene <- bind_cols(TXNAME=feat$ensembl_transcript_id_version,
                     GENEID=feat$ensembl_gene_id_version)
### Load quants
files <- list.files("~/Downloads/bcell2/", pattern = "quant.sf", full.names = T,recursive = T)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")
### Rename columns
nms <- basename(gsub("\\/quant.sf", "", files))

mat_raw <- txi$counts
colnames(mat_raw) <- gsub("_cds","", nms)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
### Rename columns
nms <- basename(gsub("\\/quant.sf", "", files))

mat <- txi$counts
colnames(mat) <- gsub("_cds","", nms)


fmat <- feat[!duplicated(feat$ensembl_gene_id_version),]
rownames(fmat) <- fmat$ensembl_gene_id_version
fmat <- fmat[rownames(mat),]
save(mat, mat_raw, fmat, file="bcell_salmon2.rda")


trsnc <- read_tsv("~/Downloads/cds/WT_CTRL1.gz_cds/quant.sf")[,1, drop = T]
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl",version = "104")
feat <- getBM(attributes=c("ensembl_transcript_id_version",
                           "ensembl_gene_id_version",
                           "external_gene_name",
                           "description"),
              filters = "ensembl_transcript_id_version",
              values = trsnc, #keys values for query
              mart = ensembl)
feat <- feat[!duplicated(feat$ensembl_transcript_id_version),]
tx2gene <- bind_cols(TXNAME=feat$ensembl_transcript_id_version,
                     GENEID=feat$ensembl_gene_id_version)
### Load quants
files <- list.files("~/Downloads/cds/", pattern = "quant.sf", full.names = T,recursive = T)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")
### Rename columns
nms <- basename(gsub("\\/quant.sf", "", files))

mat_raw_cds <- txi$counts
colnames(mat_raw_cds) <- gsub(".gz_cds","", nms)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
### Rename columns
nms <- basename(gsub("\\/quant.sf", "", files))

mat_cds <- txi$counts
colnames(mat_cds) <- gsub(".gz_cds","", nms)


fmat <- feat[!duplicated(feat$ensembl_gene_id_version),]
rownames(fmat) <- fmat$ensembl_gene_id_version
fmat <- fmat[rownames(mat_cds),]
fmat_cds <- fmat
save(mat_cds, mat_raw_cds, fmat_cds, file="bcell_salmon_cds.rda")