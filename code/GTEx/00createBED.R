library(megadepth)
library(recount3)
library(furrr)
library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Gviz)
library(annotatr)
genesTrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene)
## Add Gene symbols
myIds <- gene(genesTrack)
myIds <- AnnotationDbi::select(org.Hs.eg.db, keys=myIds, keytype="ENTREZID",columns = "SYMBOL")
#myIds <- paste(myIds[,2], myIds[,1], sep="\n")
symbol(genesTrack) <- myIds$SYMBOL
#genesTrack <-genesTrack[!is.na(symbol(genesTrack))]
apadb <- read_tsv( "~/Downloads/output_utrs_hg38.bed", col_names = F)
names(apadb) <- c("chr", "start", "end", "trash", "score", "strand", "gene_name")
apadb <- makeGRangesFromDataFrame(apadb, keep.extra.columns = T)
apadb[strand(apadb) == "+"] <- GenomicRanges::shift(flank(apadb[strand(apadb) == "+"], -50, start = F), -50)
apadb[strand(apadb) == "-"] <- GenomicRanges::shift(flank(apadb[strand(apadb) == "-"], -50, start = F), 50)
apadb$source <- "apadb2"

apadb <- sort(apadb)

apadb$annotation <- NA
#genesTrack <-genesTrack[!is.na(symbol(genesTrack))]
gr <- genesTrack@range
gr[strand(gr) == "+"] <- GenomicRanges::shift(gr[strand(gr) == "+"], -50)
gr[strand(gr) == "-"] <- GenomicRanges::shift(gr[strand(gr) == "-"], 50)
gr.cds <- gr[gr$feature == "CDS"]
gr.utr5 <- gr[gr$feature == "utr5"]
gr.utr3 <- gr[gr$feature == "utr3"]
gr.te <- gr.utr3[order(gr.utr3$exon)]
spl <- split(gr.te$exon, gr.te$transcript)

tes <- map_chr(spl, function(x){
    i <- which.max(as.numeric(gsub(".+_", "", x)))
    x[i]
})
gr.te <- gr.te[gr.te$exon %in% tes]
apadb$isTE <- FALSE

ht <- findOverlaps(apadb,gr.utr3)
apadb[seq_along(apadb) %in% unique(queryHits(ht)),]$annotation <- "3UTR"

ht <- findOverlaps(apadb,gr.cds)
apadb[seq_along(apadb) %in% unique(queryHits(ht)),]$annotation <- "CDS"

introns <- build_annotations(genome = 'hg38', annotations = c("hg38_genes_introns"))
ht <- findOverlaps(apadb,introns)
apadb[seq_along(apadb) %in% unique(queryHits(ht)),]$annotation <- "IN"


ht <- findOverlaps(apadb,gr.utr5)
apadb[seq_along(apadb) %in% unique(queryHits(ht)),]$annotation <- "5UTR"

ht <- findOverlaps(apadb,gr.te)
apadb[seq_along(apadb) %in% unique(queryHits(ht)),]$isTE <- TRUE
#exons <- build_annotations(genome = 'hg38', annotations = c("hg38_genes_exons"))
exons <- genesTrack@range

ht <- findOverlaps(apadb,exons, ignore.strand = F)
ht <- tibble(as_tibble(ht), gene_name=exons$symbol[subjectHits(ht)])
# how many genes overlaps with PA site
ht <- ht %>% group_by(queryHits) %>%
    mutate(n=length(unique(gene_name))) %>%
    mutate(drop= (n>1))
### remove double hits and then repeat with strand
keep <- unique(ht[ht$drop == FALSE,]$queryHits)
apadb <- apadb[keep]

ht <- findOverlaps(apadb,exons, ignore.strand = T)
ht <- tibble(as_tibble(ht), gene_name=exons$symbol[subjectHits(ht)])


# how many genes overlaps with PA site
ht <- ht %>% group_by(queryHits) %>%
    mutate(n=length(unique(gene_name))) %>%
    mutate(drop= (n>1))

ht1 <- ht[ht$queryHits %in% which(apadb$annotation == "3UTR") & ht$n > 1,]
ht1$strand1 <- as.character(strand(apadb[ht1$queryHits]))
ht1$strand2 <- as.character(strand(exons[ht1$subjectHits]))
new <- GRangesList(map(unique(ht1$queryHits), function(i){
    sst <- ht1[ht1$queryHits == i,]
    std <- unique(sst$strand1)
    if ( std == "+") {
        x <- flank(exons[sst[sst$strand2 == "-",]$subjectHits], start = F, 50)[1]
        strand(x) <- "+"
        x@elementMetadata <- apadb[sst[sst$strand2 == "+",]$queryHits][1]@elementMetadata
        x
    } else {
        x <- flank(exons[sst[sst$strand2 == "+",]$subjectHits], start = F, 50)[1]
        strand(x) <- "-"
        x@elementMetadata <- apadb[sst[sst$strand2 == "-",]$queryHits][1]@elementMetadata
        x
    }
}))
new <- unlist(new)
apadb[unique(ht1$queryHits)] <- new

ht <- findOverlaps(apadb,exons, ignore.strand = T)
ht <- tibble(as_tibble(ht), gene_name=exons$symbol[subjectHits(ht)])
# how many genes overlaps with PA site
ht <- ht %>% group_by(queryHits) %>%
    mutate(n=length(unique(gene_name))) %>%
    mutate(drop= (n>1 | duplicated(queryHits)))
### remove double hits and then repeat with strand
keep <- ht[ht$drop == FALSE,]$queryHits
apadb <- apadb[keep]

exons <- gr
exons <- exons[!is.na(exons$symbol)]
exons <- exons[!grepl("^MIR", exons$symbol)]

ht <- findOverlaps(apadb,exons)
ht <- tibble(as_tibble(ht), gene_name=exons$symbol[subjectHits(ht)])

# how many genes overlaps with PA site
ht <- ht %>% group_by(queryHits) %>%
    mutate(n=length(unique(gene_name))) %>%
    filter(n==1 & !duplicated(queryHits))
apadb$SYMBOL <- NA
apadb$SYMBOL[ht$queryHits] <- ht$gene_name

apadb <- apadb[!is.na(apadb$SYMBOL) & !is.na(apadb$annotation),]
apadb <- apadb[!apadb$annotation %in% c("5UTR"),]
apadb <- sort(apadb)
names(apadb) <- apadb$trash

table(sort(apadb) == apadb)

apadb$SYMBOL<- apadb$gene_name
idx <- c(names(apadb[strand(apadb) == "+"]), rev(names(apadb[strand(apadb) == "-"])))
apadb <- apadb[idx,]
ids <- tibble(symbol=apadb$SYMBOL) %>% group_by(symbol) %>%
    mutate(ids=paste(symbol, sprintf("%02d", 1:n()), sep="_"))
ids <- ids$ids
apadb$ids <- NA
apadb$ids <- ids
names(apadb) <- ids
apadb <- apadb[sort(names(apadb)),]
write_tsv(data.frame(apadb)[,c("seqnames", "start", "end", "ids", "score", "strand", "SYMBOL", "annotation", "isTE", "qapa_id")], file = "data/QAPA_up_hg38.bed", col_names = F)
