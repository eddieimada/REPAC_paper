library(tidyverse)
library(edgeR)
# Load QAPA results
qapa <- read_csv("../APA/text/paper/QAPA_results.csv")
#load gene level counts and filter low expressed genes
load("objs/bcell_salmon.rda")
dge <- DGEList(mat, genes = fmat, group = se$sra.sample_title)
keep <- filterByExpr(dge, group = se$sra.sample_title, min.count= 30, min.total.count=10)
table(keep)
gns <- gsub("\\..+", "", fmat[keep,]$ensembl_gene_id_version)
qapa <- qapa[qapa$Gene %in% gns,]

# remove UTRs with very low expression
qapaX <-qapa %>% dplyr::select(ends_with("TPM"))
keep <- rowSums(qapaX[,1:4] >= 10) == 4 | rowSums(qapaX[,5:8] >= 10) == 4 
table(keep)
qapa<- qapa[keep,]

# compute delta PAU
qapaA <-qapa %>% dplyr::select(contains("CTRL") & ends_with("PAU"))
qapaA <- rowMeans(qapaA,na.rm = T)
qapaB <-qapa %>% dplyr::select(contains("LPS") & ends_with("PAU"))
qapaB <- rowMeans(qapaB, na.rm = T)
qapa_score <- qapaA-qapaB

qapa$pau_diff <- qapa_score
qapa.p <- qapa[grepl("P$", qapa$APA_ID),]
write_csv(qapa.p, file="qapa_ppaus.csv")
qapa.final <- qapa[qapa_score <= -20 | qapa_score >= 20,]
qapa.final <- qapa.final[order(abs(qapa.final$pau_diff), decreasing = T),]
# select PPAUs
qapa.final.p <- qapa.final[grepl("P$", qapa.final$APA_ID),]

write_csv(qapa.final.p, file="text/paper/qapa_final.csv")
