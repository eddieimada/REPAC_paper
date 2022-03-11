library(GenomicRanges)
library(GenomicFeatures)
library(data.table)
library(tidyverse)

# This script creates a new GTF with 5000 randomly selected genes with a long and a short isoform.
anno <- fread("~/Dropbox (MechPred)/Databases/GTFs/hg38/GENCODE/gencode.v37.annotation.gtf.gz")
names(anno)[1] <- "V1"
gtf <- anno
#Split annotations
tmp <- gtf$V9
tmp <- gsub(" level [0-9];", "", tmp)
sp <- strsplit(tmp, "; ")
tmp <- map(sp, ~{
    lgl <- grepl("ont |tag ", .x)
    tmp <- .x[!lgl]
    val <- gsub(".+ \"(.+)\"", "\\1", tmp)
    id <- gsub("(.+) .+", "\\1", tmp)
    set_names(val,id)
}
)
tmp2<- bind_rows(tmp)
gtf <- gtf[,-9]
gtf <- cbind(gtf, tmp2)
names(gtf)[1] <- "V1"

# select protein coding genes
gtf <- gtf[gtf$V1 %in% paste0("chr", c(1:22, "X", "Y")),]
gtf <- gtf[gtf$gene_type == "protein_coding",]
gtf <- gtf[gtf$transcript_type == "protein_coding" | is.na(gtf$transcript_type),]

# Select the longest transcript
tmp <- split(gtf, gtf$gene_id)
tmp <- map(tmp, ~{
    .x %>% 
        group_by(transcript_id) %>%
        filter(!is.na(transcript_id)) %>% 
        transmute(len=sum(V5[V3 == "exon"] - V4[V3 == "exon"])) %>%
        ungroup() %>%
        transmute(id=transcript_id[which.max(len)])
})
sel <- map_chr(tmp,~ ifelse(nrow(.x) != 0, unique(.x$id), NA))
sel <- sel[!is.na(sel)]
gtf <- gtf[gtf$transcript_id %in% sel | is.na(gtf$transcript_type),]

# Keep only genes 5000bp away from each other
colnames(gtf)[c(1,4,5,7)] <- c("chr", "start", "stop", "strand")
gtf.gr <- makeGRangesFromDataFrame(gtf[gtf$V3=="gene",], keep.extra.columns = T)
start(gtf.gr) <- case_when(as.character(strand(gtf.gr)) == "-" ~ start(gtf.gr) - 5000L,
                           TRUE ~ start(gtf.gr))
end(gtf.gr) <- case_when(as.character(strand(gtf.gr)) == "+" ~ end(gtf.gr) + 5000L,
                         TRUE ~ end(gtf.gr))
gns.u <- gtf.gr[countOverlaps(gtf.gr,gtf.gr) == 1]$gene_id
gtf <- split(gtf, gtf$gene_id)
gtf <- gtf[names(gtf) %in% gns.u]
set.seed(666)
idx <- sample(seq_along(gtf), size = 5000, replace = F)
gtf <- gtf[idx]
gtf <- bind_rows(gtf)
colnames(gtf)[c(1,4,5,7)] <- c("V1", "V4", "V5", "V7")

ids.n <- paste(gtf$V1, gtf$V2, gtf$V3, gtf$V4, gtf$V5,gtf$V7, gtf$gene_id, gtf$transcript_id, sep = "_")
tmp <- anno[,-9]
tmp <- cbind(tmp, tmp2)
names(tmp)[1] <- "V1"
ids <- paste(anno$V1, anno$V2, anno$V3, anno$V4, anno$V5,anno$V7, tmp$gene_id, tmp$transcript_id, sep = "_")
anno <- cbind(anno,gene_id=tmp2$gene_id)
anno.sel <- anno[ids %in% ids.n,]

saveRDS(anno.sel, file="anno_sel.rds")
# Expand last exon by 1000bp
anno.exp1 <- anno.sel %>%
    group_by(gene_id) %>%
    
    mutate(V5 = case_when(V7 == "+" & V3 == "exon" & 
                              V5 == max(V5[V3 == "exon"]) ~ V5 + 1000L,
                          TRUE ~ V5),
           V4 = case_when(V7 == "-" & V3 == "exon" & 
                              V4 == min(V4[V3 == "exon"])~ V4 - 1000L,
                          TRUE ~ V4)) %>%
    mutate(V5 = case_when(V7 == "+" & V3 %in% c("gene", "transcript") ~ max(V5[V3 == "exon"]),
                          TRUE ~ V5),
           V4 = case_when(V7 == "-" & V3 %in% c("gene", "transcript") ~ min(V4[V3 == "exon"]),
                          TRUE ~ V4))
gns <- unique(anno.sel$gene_id)
anno.sel$V9 <- gsub("ENST", "S",anno.sel$V9)
anno.exp1$V9 <- gsub("ENST", "L",anno.exp1$V9)
anno.sel <- anno.sel[anno.sel$V3 %in% c("gene", "transcript", "exon", "start_codon", "CDS", "stop_codon"),]
anno.exp1 <- anno.exp1[anno.exp1$V3 %in% c("gene", "transcript", "exon", "start_codon", "CDS", "stop_codon"),]

# create new gtf with both isoforms
tmp1 <- anno.sel[anno.sel$V3 != "gene",]
tmp <- bind_rows(anno.exp1, tmp1)
tmp$V9 <- factor(tmp$V9, levels = unique(anno.exp1$V9))
tmp$V1 <- factor(tmp$V1, levels = c(paste0("chr", c(1:22, "X", "Y"))), ordered = T)
tmp$gene_id <- factor(tmp$gene_id, levels = unique(tmp$gene_id), ordered = T)
tmp$V10 <- gsub("gene_type .+", "", tmp$V9)
tmp <- tmp %>%
    arrange(V9)
setC<- tmp[,c(-10,-11)]

write.table(setC, file = "setC_2iso_5k1k_new.gtf",col.names = F, row.names = F, quote = F, sep = "\t")

# create gtf with 50bp windows from PA site
pa <- tmp %>%
    group_by(gene_id) %>%
    filter(V3 == "gene" | V3 == "transcript" |
               (V3 == "exon" & V7 == "+" & V5 == V5[V3 == "gene"]) |
               (V3 == "exon" & V7 == "-" & V4 == V4[V3 == "gene"])) %>%
    mutate(V4 = case_when(V7 == "+" & V3 %in% c("exon", "transcript") ~ V5 - 100L,
                          TRUE ~ V4),
           V5 = case_when(V7 == "+" & V3 %in% c("exon", "transcript") ~ V5 - 50L,
                          TRUE ~ V5),
           V5 = case_when(V7 == "-" & V3 %in% c("exon", "transcript") ~ V4 + 100L,
                          TRUE ~ V5),
           V4 = case_when(V7 == "-" & V3 %in% c("exon", "transcript") ~ V4 + 50L,
                          TRUE ~ V4)
           
    ) 

pa<- pa[,c(-10,-11)]
write.table(pa, file = "pa_sim_2iso_5k1k_new.gtf",col.names = F, row.names = F, quote = F, sep = "\t")

