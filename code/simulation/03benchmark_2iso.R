library(caret)
library(tidyverse)
library(ROCit)
library(ggplot2)

# This scripts documents the benchmarking process. The results from QAPA and DaPars are included in this repository. For details on how they were generated, please refer to the publication.

# load profile and obtain true values
proA <- read_tsv("parA.pro", col_names = F)
idx <- c(rep("1", 1250), rep("2", 1250), rep("3", 2500))
set.seed(777)
idx <- sample(idx)
pct1 <- unlist(map(idx, ~{
    if (.x == 1) {
        r <- sample(seq(0.1,0.35, by=0.01), 1)
        c(1-r,r)
    } else if (.x==2) {
        r <- sample(seq(0.1,0.35, by=0.01), 1)
        c(r,1-r)
    } else {
        c(0.5,0.5)
    }
}))

truth <- cbind.data.frame(id=proA$X2, pct=pct1)
tt <- data.frame(truth) %>%
    mutate(tt=case_when(pct == 0.5 ~ "NORMAL",
                        pct > 0.35 ~ "LONG",
                        pct <= 0.35 ~ "SHORT"))
tt <- set_names(tt$tt, tt$id) 

# read REPAC results
simRes <- read_csv("simStats2iso_weights_test.csv")
gp <- simRes 
gp$pred <- case_when(simRes$cFC >= 0.2 & simRes$adj.p.val <= 0.1 ~ "SHORT",
              simRes$cFC <= -0.2 & simRes$adj.p.val <= 0.1 ~ "LONG",
              TRUE ~"NORMAL")
tt2 <- tt[grepl("^L", names(tt))]
names(tt2) <-  gsub(".(.+)", "\\1_02", names(tt2))
res <- cbind(gp, tt=tt2[gp$ID])
table(res$pred == res$tt)
pred.res <- factor(res$pred, levels= c("LONG","SHORT", "NORMAL"))
trt <- factor(res$tt, levels= c("LONG","SHORT", "NORMAL"))
caret::confusionMatrix(pred.res, trt)
class <- ifelse(trt == "NORMAL", 0, 1)
pred <- ifelse(pred.res == "NORMAL", 0, 1)
# ROC curve
rocit_corepad <- rocit(score = abs(gp$cFC), class = class)
summary(rocit_corepad)
################################################################
# read QAPA results
pau <- read_tsv("~/simulation/QAPA_final.tsv")
pau.proximal <- pau %>% filter(grepl("_P$", APA_ID)) 
qapaA <-data.frame(pau.proximal) %>% dplyr::select(starts_with("simA") & ends_with("PAU"))
qapaA <- rowMeans(qapaA)
qapaB <-data.frame(pau.proximal) %>% dplyr::select(starts_with("simB") & ends_with("PAU"))
qapaB <- rowMeans(qapaB)
qapa_score <- qapaA-qapaB
qapa <- qapaA-qapaB

qapa <- case_when(qapa <= -20 ~"LONG",
                  qapa >= 20 ~ "SHORT",
                  TRUE~"NORMAL")

qapa <- cbind.data.frame(pau.proximal$Transcript, gp=qapa)
qapa.res <- cbind.data.frame(gp=qapa$gp, tt=tt[qapa$`pau.proximal$Transcript`])

table(qapa.res$gp == qapa.res$tt)

pred.res <- factor(qapa.res$gp, levels= c("LONG","SHORT", "NORMAL"))
trt <- factor(qapa.res$tt, levels= c("LONG","SHORT", "NORMAL"))
caret::confusionMatrix(pred.res, trt)
class <- ifelse(trt == "NORMAL", 0, 1)
rocit_qapa <- rocit(score = abs(qapa_score), class = class)
summary(rocit_qapa)
#################################################
# read Dapars results
dp <-read_tsv(file = "DaPars_Result_All_Prediction_Results.txt", col_names = T)
dp$ID <- map_chr(strsplit(dp$Gene, "\\|"),~.x[[1]])
cutoff <- 0.2
dapars <- dp %>%
    mutate(SYMBOL=gsub("S|L", "", ID)) %>%
    group_by(SYMBOL) %>%
    filter(n() ==2 ) %>%
    ungroup()%>%
    mutate(gp=case_when(adjusted.P_val > 0.1 |
                            PDUI_Group_diff > cutoff*-1 & PDUI_Group_diff < cutoff ~ "NORMAL",
                        adjusted.P_val <= 0.1 & PDUI_Group_diff >= cutoff ~ "LONG",
                        adjusted.P_val <= 0.1 & PDUI_Group_diff <= cutoff*-1 ~ "SHORT")) %>%
    filter(grepl("^L", ID))

dap <- cbind(dapars, tt=tt[dapars$ID])

pred.dap <- factor(dap$gp, levels= c("LONG","SHORT", "NORMAL"))
trt <- factor(dap$tt, levels= c("LONG","SHORT", "NORMAL"))
caret::confusionMatrix(pred.dap, trt)
aa <- a[grepl("^L", a$ID),]
class <- ifelse(tt[aa$ID] == "NORMAL", 0, 1)
rocit_dapars <- rocit(score = abs(aa$PDUI_Group_diff), class = class)
summary(rocit_dapars)
 
#################################################
# Plot ROC curves
corepad <- cbind.data.frame(TPR=rocit_corepad$TPR,FPR=rocit_corepad$FPR)
qapa <- cbind.data.frame(TPR=rocit_qapa$TPR,FPR=rocit_qapa$FPR)
dapars <- cbind.data.frame(TPR=rocit_dapars$TPR,FPR=rocit_dapars$FPR)
l <- list(corepad,qapa,dapars)
names(l) <- c("REPAC (AUC=0.9989)", "QAPA   (AUC=0.9967)", "DaPars (AUC=0.8820)")
df <- bind_rows(l, .id="Method")

p <- ggplot(df, aes(x=FPR, y=TPR, color=Method)) +
    geom_line(lwd = 1.1, alpha=0.7) +
    geom_abline(lty =2, lwd =1.2) +
    ylab("Sensitivity (TPR)") +
    xlab("1 - Specificity (FPR)") +
    scale_x_continuous(expand = c(0.01, 0.03), limits = c(0, 1)) +
    scale_y_continuous(expand = c(0, 0.01), limits = c(0, 1)) +
    theme_bw() +
    theme(text = element_text(size=14), 
          legend.text = element_text(size=14),
          axis.text = element_text(size=14), legend.title = element_blank(),
          legend.justification=c(1,0), legend.position=c(1,0),
          legend.background = element_rect(fill="transparent")) 

ggsave(p, file="figs/AUC_simulation.png", width = 6, height = 6, dpi = 330)
