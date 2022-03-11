library(derfinder)
library(SummarizedExperiment)
library(edgeR)
library(data.table)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(org.Mm.eg.db)
library(Gviz)
library(tidyverse)
library(annotatr)


se <- read_rds("objs/SRP048707_APA_qapa.rds")
se <- se[,grepl("WT", se$sra.sample_title)]

gr <- rowRanges(se)
gr$ids <- names(gr)
gr$gene_name  <- gr$SYMBOL

# load gene tracks
genesTrack <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, 
                              name="", col = "black",
                              fill="orange", col.line="black", lwd=0.8, size =1, stacking = "dense")
genesTrack <- genesTrack[genesTrack@range$feature != "ncRNA"]
## Add Gene symbols
myIds <- gene(genesTrack)
myIds <- AnnotationDbi::select(org.Mm.eg.db, keys=myIds, keytype="ENTREZID",columns = "gene_name")
symbol(genesTrack) <- myIds$gene_name
genesTrack <-genesTrack[!is.na(symbol(genesTrack))]

exons <- genesTrack@range
exons <- exons[!is.na(exons$symbol)]
exons <- exons[!grepl("^MIR", exons$symbol)]

#load results
tb.g.corr <- read_csv("bcell.csv")

tb.g.corr <- tb.g.corr %>%
    arrange(abs(p.val)) %>%
    filter(!duplicated(Ref))

tb.ord <- tb.g.corr %>%
    filter(abs(cFC) >= 0.25 & adj.p.val <= 0.05 ) %>%
    arrange(adj.p.val)

# select top 10 DPU
gns <- tb.ord
gns <- gns[!duplicated(gns$gene_name),]
gns1 <- head(gns$ID,10)
gns2 <- head(gns$Ref,10)
bws <- se$BigWigURL
names(bws) <- colnames(se)

# Retrieve coverage for select samples and ranges
exons <- exons[seqnames(exons) %in% c(paste0("chr", c(1:22, "X", "Y", "M")))]
tt <- sapply(seq_along(gns1), function(nm){
    gn <- gr[names(gr) %in% c(gns1[nm], gns2[nm])]
    trnsc <- GenomicRanges::reduce(exons[exons$symbol %in% 
                                             gsub("_.+", "", c(gns1[nm], gns2[nm]))])
    trnsc <- trnsc[trnsc %over% gn]
})

tt <- map(tt, function(x){
    GRanges(seqnames = seqnames(x)[1],
            ranges=IRanges(
                start = min(start(x)),
                end = max(end(x))),
            strand=strand(x)[1])
})
names(tt) <- gns1

regionCov <- getRegionCoverage(regions = unlist(GRangesList(tt)), files = bws,
                               verbose = T)
regionCov <- map(regionCov, function(y) {
    colnames(y) <- names(bws)
    y})
object.size(regionCov)

names(regionCov) <- gns1
# rescale for plotting
regionCov2 <- map(regionCov, function(x){
    x <- as.matrix(x)
    sweep(x, 2, colMaxs(x), "/")
})

sapply(seq_along(names(regionCov2)), function(i){
    nm <- names(regionCov2)[[i]]
    png(paste0("figs/mm/new_t_bcell_",sprintf("%02d", i),"_", nm,".png"), width=2000, height=1000, res = 330)
    print(nm)
    gn <- gr[nm]
    chr <- as.character(seqnames(gn)[1])
    gn <- gr[names(gr)  %in% c(gns1[i], gns2[i])]
    trnsc <- GenomicRanges::reduce(exons[exons$symbol %in% 
                                             gsub("_.+", "", c(gns1[i], gns2[i]))])
    trnsc <- trnsc[trnsc %over% gn]
    wdt <- max(end(trnsc)) - min(start(trnsc))
    gene.model <- GRanges(seqnames = seqnames(trnsc)[1],
                          ranges=IRanges(
                              start = min(start(trnsc)),
                              end = max(end(trnsc))),
                          strand=strand(trnsc)[1]
    )
    keep <- match(colnames(regionCov2[[i]]),
                  colnames(se))
    pheno <- colData(se)[keep,]
    groups <- factor(pheno$sra.sample_title)
    groups <- relevel(groups, ref = "WT_Ex vivo")
    mat <- as.matrix(regionCov2[[i]])
    
    mat.cpm <- as.data.frame(mat)
    chr <- as.character(seqnames(gn)[1])
    
    st <-(start(gene.model) + 0:(nrow(mat.cpm)-1))
    fst <- st[1]
    lst <- st[length(st)]
    chr <- as.character(seqnames(gn)[1])
    ann <- data.frame(chr=rep(chr,length(st)),
                      start=st,
                      end=st
    )
    ann <- bind_cols(ann,mat.cpm)
    apa.gr <- makeGRangesFromDataFrame(ann, keep.extra.columns = T)
    ref <- rownames(pheno)[pheno$sra.sample_title== "WT_Ex vivo"]
    gs.ref <- factor(ref)
    ref <- apa.gr[,ref]
    contr <- rownames(pheno)[pheno$sra.sample_title == "WT_LPS"]
    gs.contr <- factor(contr)
    contr <- apa.gr[,contr]
    ### Create tracks
    contr.max <- 1
    ref.max <- 1

    ref2 <-  DataTrack(apa.gr[,pheno$sra.sample_title == "WT_Ex vivo"], baseline = 0, col.baseline="black", 
                       genome="hg38", type = c("h"),  col =c("red"), 
                       alpha =1,
                       showSampleNames=T, name=gsub("_.+", "", nm),size=2, cex.sampleNames = 1)
    ref3 <-  DataTrack(apa.gr[,pheno$sra.sample_title == "WT_LPS"], baseline = 0, col.baseline="black", 
                       genome="hg38", type = c("h"),  col =c("black"), 
                       alpha =1,
                       showSampleNames=T, name=gsub("_.+", "", nm),size=2, cex.sampleNames = 1)
    myChr <- as.character(seqnames(gn))[1]
    chromosome(genesTrack) <- as.character(myChr)

    gt <- genesTrack[as.character(strand(genesTrack@range)) == as.character(strand(trnsc))[1]]
    plotTracks(c(gt, ref2),
               add=FALSE, from=fst - (wdt * 0.05), aggregateGroups = TRUE, aggregation = "mean", sizes = c(0.5,2),
               to=lst + (wdt * 0.05), chromosome=myChr, grid=F,cex.axis = 1,
               stacking = "squish",
               background.title = "white", col.axis ="black", fontcolor.title = "black"
    )
    plotTracks(c(gt, ref3),
               add=T, from=fst - (wdt * 0.05), aggregateGroups = TRUE, aggregation = "mean", sizes = c(2,2),
               to=lst + (wdt * 0.05), chromosome=myChr, grid=F,cex.axis = 1,
               stacking = "squish",
               background.title = "white", col.axis ="black", fontcolor.title = "black"
    )

    dev.off()
})
