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
truth <- read_csv("interaction_truth.csv")
tt <- data.frame(truth) %>%
    mutate(tt=case_when(idx == 1 ~ "SHORT",
                        idx == 2 ~ "LONG",
                        idx == 3 ~ "LONG",
                        idx == 4 ~ "SHORT",
                        idx >= 5 ~ "NORMAL"))
tt <- set_names(tt$tt, gsub("ENSG", "",tt$gns))

# read REPAC results

gp <- read_csv("interaction_sim.csv")
gp$ID <- gsub("_.+", "", gp$Site)
gp$pred <- case_when(gp$cFC <= -0.2 & gp$P.Value <= 0.1 ~ "SHORT",
              gp$cFC >= 0.2 & gp$P.Value <= 0.1 ~ "LONG",
              TRUE ~"NORMAL")
res <- cbind(gp, tt=tt[gp$ID])
table(res$pred == res$tt)
pred.res <- factor(res$pred, levels= c("LONG","SHORT", "NORMAL"))
trt <- factor(res$tt, levels= c("LONG","SHORT", "NORMAL"))
c.int <- caret::confusionMatrix(pred.res, trt)
cm.int <- t(c.int$byClass)
class <- ifelse(trt == "NORMAL", 0, 1)
# ROC curve
rocit_int <- rocit(score = abs(gp$cFC), class = class)
summary(rocit_int)
################################################################
truth <- read_csv("2iso_truth.csv")
tt <- data.frame(truth) %>%
    mutate(tt=case_when(idx == 1 ~ "SHORT",
                        idx == 2 ~ "LONG",
                        idx == 3 ~ "NORMAL"))
tt <- set_names(tt$tt, gsub("ENSG", "",tt$gns))

# read REPAC results
gp <- read_csv("50bp.csv")
gp$ID <- gsub("_.+", "", gp$Site)
gp$pred <- case_when(gp$cFC <= -0.2 & gp$P.Value <= 0.1 ~ "SHORT",
                     gp$cFC >= 0.2 & gp$P.Value <= 0.1 ~ "LONG",
                     TRUE ~"NORMAL")
res <- cbind(gp, tt=tt[gp$ID])
table(res$pred == res$tt)
pred.res <- factor(res$pred, levels= c("LONG","SHORT", "NORMAL"))
trt <- factor(res$tt, levels= c("LONG","SHORT", "NORMAL"))
c.rz <- caret::confusionMatrix(pred.res, trt)
cm.rz <- t(c.rz$byClass)
class <- ifelse(trt == "NORMAL", 0, 1)
# ROC curve
rocit_rz <- rocit(score = abs(gp$cFC), class = class)
summary(rocit_rz)
 
#################################################
truth <- read_csv("2iso_truth.csv")
tt <- data.frame(truth) %>%
    mutate(tt=case_when(idx == 1 ~ "SHORT",
                        idx == 2 ~ "LONG",
                        idx == 3 ~ "NORMAL"))
tt <- set_names(tt$tt, gsub("ENSG", "",tt$gns))

# read REPAC results
gp <- read_csv("polyA.csv")
gp$ID <- gsub("_.+", "", gp$Site)
gp$pred <- case_when(gp$cFC <= -0.2 & gp$P.Value <= 0.1 ~ "SHORT",
                     gp$cFC >= 0.2 & gp$P.Value <= 0.1 ~ "LONG",
                     TRUE ~"NORMAL")
res <- cbind(gp, tt=tt[gp$ID])
table(res$pred == res$tt)
pred.res <- factor(res$pred, levels= c("LONG","SHORT", "NORMAL"))
trt <- factor(res$tt, levels= c("LONG","SHORT", "NORMAL"))
c.pa <- caret::confusionMatrix(pred.res, trt)
cm.pa <- t(c.pa$byClass)
class <- ifelse(trt == "NORMAL", 0, 1)
# ROC curve
rocit_pa <- rocit(score = abs(gp$cFC), class = class)
summary(rocit_pa)

#################################################
truth <- read_csv("2iso_truth.csv")
tt <- data.frame(truth) %>%
    mutate(tt=case_when(idx == 1 ~ "SHORT",
                        idx == 2 ~ "LONG",
                        idx == 3 ~ "NORMAL"))
tt <- set_names(tt$tt, gsub("ENSG", "",tt$gns))

# read REPAC results
gp <- read_csv("singleEnd.csv")
gp$ID <- gsub("_.+", "", gp$Site)
gp$pred <- case_when(gp$cFC <= -0.2 & gp$P.Value <= 0.1 ~ "SHORT",
                     gp$cFC >= 0.2 & gp$P.Value <= 0.1 ~ "LONG",
                     TRUE ~"NORMAL")
res <- cbind(gp, tt=tt[gp$ID])
table(res$pred == res$tt)
pred.res <- factor(res$pred, levels= c("LONG","SHORT", "NORMAL"))
trt <- factor(res$tt, levels= c("LONG","SHORT", "NORMAL"))
c.pa <- caret::confusionMatrix(pred.res, trt)
cm.pa <- t(c.pa$byClass)
class <- ifelse(trt == "NORMAL", 0, 1)
# ROC curve
rocit_pa <- rocit(score = abs(gp$cFC), class = class)
summary(rocit_pa)

