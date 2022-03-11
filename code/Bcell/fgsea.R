library(tidyverse)
library(fgsea)
library(BiocParallel)
library(ggplot2)
library(RColorBrewer)
library(ggnewscale)
library(ComplexHeatmap)
library(compositions)
register(SerialParam())

se <- readRDS("objs/SRP048707_APA_qapa.rds")
tb.g.corr <- read_csv("bcell.csv")
top.t <- tb.g.corr %>%
    group_by(gene_name) %>%
    filter(1:n() == which.max(abs(cFC)))

rnk.apa <- setNames(top.t$t,toupper(top.t$gene_name))

qapa.p <- read_csv("qapa_ppaus.csv")
rnk.qapa <- setNames(qapa.p$pau_diff,toupper(qapa.p$Gene_Name))


set.seed(777)
BP <- gmtPathways("~/Dropbox (MechPred)/Databases/MSigDB/v7.4/c5.go.bp.v7.4.symbols.gmt.txt")

enrich.apa <- fgseaMultilevel(BP, rnk.apa, nproc = 8, minSize = 30, maxSize = 80, eps=0, nPermSimple = 10000, scoreType = "neg") 

apa.sig <- enrich.apa[enrich.apa$pval <= 0.05,]
slim1 <- collapsePathways(apa.sig[order(apa.sig$pval),], BP[names(BP) %in% apa.sig$pathway], rnk.apa, nperm=1000)
apa.sig <- apa.sig %>% filter(pathway %in% slim1$mainPathways)
apa.sig$padj <- p.adjust(apa.sig$pval, method = "BH")

apa.sig <- apa.sig[order(apa.sig$pval),]

write_csv(apa.sig, file="text/paper/GOBP_REPAC.csv")

enrich.qapa <- fgseaMultilevel(BP, rnk.qapa, nproc = 8, minSize = 30, maxSize = 80, eps=0, nPermSimple = 10000, scoreType = "neg") 

qapa.sig <- enrich.qapa[enrich.qapa$pval <= 0.05,]
slim <- collapsePathways(qapa.sig[order(qapa.sig$pval),], BP[names(BP) %in% qapa.sig$pathway], rnk.qapa, nperm=1000)
qapa.sig <- qapa.sig %>% filter(pathway %in% slim$mainPathways)
qapa.sig$padj <- p.adjust(qapa.sig$pval, method = "BH")

qapa.sig <- qapa.sig[order(qapa.sig$padj),]

write_csv(qapa.sig, file="text/paper/GOBP_QAPA.csv")



secretion  <- unlist(enrich.apa[enrich.apa$pathway == "GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",]$leadingEdge)
interferon  <- unlist(enrich.apa[enrich.apa$pathway == "GOBP_RESPONSE_TO_TYPE_I_INTERFERON",]$leadingEdge)

pa.counts <- assays(se)[[1]]

sec <- tb.g.corr[match(secretion, toupper(tb.g.corr$gene_name)),] %>% 
    filter(p.val <= 0.05) %>%
    arrange(cFC)

sec.df <- map2(sec$Ref, sec$ID, function(ref,id){
    as.numeric(ilr(t(pa.counts[c(ref,id),])))
})
names(sec.df) <- sec$gene_name


int <- tb.g.corr[match(interferon, toupper(tb.g.corr$gene_name)),] %>% 
    filter(p.val <= 0.05) %>%
    arrange(cFC)

int.df <- map2(int$Ref, int$ID, function(ref,id){
    as.numeric(ilr(t(pa.counts[c(ref,id),])))
})
names(int.df) <- int$gene_name

sec.df <- bind_rows(sec.df)
int.df <- bind_rows(int.df)
comps <- t(bind_cols(sec.df,int.df)
)
comps <- t(scale(t(comps)))

load("bcell_salmon.rda")
dge <- DGEList(mat, genes = fmat, group = se$sra.sample_title)
keep <- filterByExpr(dge, group = se$sra.sample_title, min.count= 5, min.total.count=10)
table(keep)
dge <- dge[keep,]
dge <- calcNormFactors(dge)
mat.cpm <- cpm(dge)
rownames(mat.cpm) <- dge$genes$external_gene_name
exps <- t(scale(t(mat.cpm[c(sec$gene_name, int$gene_name),])))

ha = HeatmapAnnotation(
    foo = anno_block(gp = gpar(fill = c(4,6)), labels = c("ACTIVATED","NAIVE"))
    )
split_h = rep(1:2, each = 4)

ra = rowAnnotation(
    foo = anno_block(gp = gpar(fill = c(2,3)), labels = c("IRE1 Response", "IFN-I Response")))
split = c(rep(1, nrow(sec)), rep(2, nrow(int)))

h1 <- Heatmap(comps, left_annotation = ra, row_split = split,
              top_annotation = ha, column_split = split_h,
        cluster_columns = T, column_title = "PPAS:DPAS Usage",
        name="h1",
        #row_labels = gsub("\\..+", "", rownames(comps)),
        show_row_names = F,
        show_column_names = F,
        heatmap_legend_param = list(title = ""),
        cluster_rows = F,
        rect_gp= gpar(col= "white", lwd=2.5),
        show_row_dend = F,
        row_names_side = "left"
        )

h2 <- Heatmap(exps, row_split = split,
              top_annotation = ha, column_split = split_h,
              cluster_columns = T,column_title = "Gene Expression",
              name="h2",
              #row_labels = gsub("\\..+", "", rownames(comps)),
              show_row_names = F,
              show_column_names = F,
              show_heatmap_legend = F,
              cluster_rows = F,
              rect_gp= gpar(col= "white", lwd=2.5),
              show_row_dend = F,
              row_names_side = "left"
)

hm <- h1+h2

png(file= "figs/mm/heatmap2.png", height = 1800, width = 2000, res = 330)
hm
dev.off()




