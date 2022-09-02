library(tidyverse)

# This scripts creates the base expression profiles for the 2 conditions simulated with simFlux. The base profiles are provided in this repository.

# Load expression profiles for the simulation
proA <- read_tsv("~/Downloads/flux-simulator-1.2.1/ribozero/parA.pro", col_names = F)
gns <- unique(gsub("S|L", "ENSG", proA$X2))
set.seed(777)
# Randomly assign a ratio of long:short isoform between 10:90 to 35:65

idx <- c(rep("1", 625), 
         rep("2", 625), 
         rep("3", 625), 
         rep("4", 625),
         rep("5", 625), 
         rep("6", 625),
         rep("7", 1250))
idx <- sample(idx)
truth <- cbind(gns,idx)
write_csv(data.frame(truth), "interaction_truth.csv")
pct1 <- map(idx, ~{
    if (.x == 1) {
        r <- sample(seq(0.1,0.15, by=0.01), 1)
        y <- sample(seq(0.3,0.35, by=0.01), 1)
        c(1-r,r,1-y,y)
    } else if (.x == 2) {
        r <- sample(seq(0.1,0.15, by=0.01), 1)
        y <- sample(seq(0.3,0.35, by=0.01), 1)
        c(r,1-r,y,1-y)
    } else if (.x == 3) {
        r <- sample(seq(0.1,0.15, by=0.01), 1)
        y <- sample(seq(0.3,0.35, by=0.01), 1)
        c(1-y,y, 1-r,r)
    } else if (.x == 4) {
        r <- sample(seq(0.1,0.15, by=0.01), 1)
        y <- sample(seq(0.3,0.35, by=0.01), 1)
        c(y,1-y,r,1-r)
    } else if (.x == 5) {
        r <- sample(seq(0.1,0.35, by=0.01), 1)
        c(1-r,r,1-r,r)
    } else if (.x == 6) {
        r <- sample(seq(0.1,0.35, by=0.01), 1)
        c(r,1-r,r,1-r)
    } else {
        c(0.5,0.5,0.5,0.5)
    }
})

p1 <- unlist(map(pct1, ~ .x[1:2]))
p2 <- unlist(map(pct1, ~ .x[3:4]))
proC <- bind_cols(proA,pct=p2)
proA <- bind_cols(proA,pct=p1)

proA <- proA %>%
    mutate(id=gsub("L|S", "", X2)) %>%
    group_by(id) %>%
    mutate(X6=round(sum(X6)*pct)) %>%
    ungroup() %>%
    mutate(X5=X6/sum(X6)) %>%
    dplyr::select(-pct,-id)

proC <- proC %>%
    mutate(id=gsub("L|S", "", X2)) %>%
    group_by(id) %>%
    mutate(X6=round(sum(X6)*pct)) %>%
    ungroup() %>%
    mutate(X5=X6/sum(X6)) %>%
    dplyr::select(-pct,-id)
# Same for the second condition
proB <- read_tsv("~/Downloads/flux-simulator-1.2.1/ribozero/parB.pro", col_names = F)

pct1 <- map(idx, ~{
    if (.x == 2) {
        r <- sample(seq(0.1,0.15, by=0.01), 1)
        y <- sample(seq(0.3,0.35, by=0.01), 1)
        c(1-r,r,1-y,y)
    } else if (.x == 1) {
        r <- sample(seq(0.1,0.15, by=0.01), 1)
        y <- sample(seq(0.3,0.35, by=0.01), 1)
        c(r,1-r,y,1-y)
    } else if (.x == 4) {
        r <- sample(seq(0.1,0.15, by=0.01), 1)
        y <- sample(seq(0.3,0.35, by=0.01), 1)
        c(1-y,y, 1-r,r)
    } else if (.x == 3) {
        r <- sample(seq(0.1,0.15, by=0.01), 1)
        y <- sample(seq(0.3,0.35, by=0.01), 1)
        c(y,1-y,r,1-r)
    } else if (.x == 6) {
        r <- sample(seq(0.1,0.35, by=0.01), 1)
        c(1-r,r,1-r,r)
    } else if (.x == 5) {
        r <- sample(seq(0.1,0.35, by=0.01), 1)
        c(r,1-r,r,1-r)
    } else {
        c(0.5,0.5,0.5,0.5)
    }
})

p1 <- unlist(map(pct1, ~ .x[1:2]))
p2 <- unlist(map(pct1, ~ .x[3:4]))
proD <- bind_cols(proB,pct=p2)
proB <- bind_cols(proB,pct=p1)

proB <- proB %>%
    mutate(id=gsub("L|S", "", X2)) %>%
    group_by(id) %>%
    mutate(X6=round(sum(X6)*pct)) %>%
    ungroup() %>%
    mutate(X5=X6/sum(X6)) %>%
    dplyr::select(-pct,-id)

proD <- proD %>%
    mutate(id=gsub("L|S", "", X2)) %>%
    group_by(id) %>%
    mutate(X6=round(sum(X6)*pct)) %>%
    ungroup() %>%
    mutate(X5=X6/sum(X6)) %>%
    dplyr::select(-pct,-id)
# Save profiles to use with simFlux
write.table(proA, file = "simA_int.pro",col.names = F, row.names = F, quote = F, sep = "\t")
write.table(proB, file = "simB_int.pro",col.names = F, row.names = F, quote = F, sep = "\t")
write.table(proC, file = "simC_int.pro",col.names = F, row.names = F, quote = F, sep = "\t")
write.table(proD, file = "simD_int.pro",col.names = F, row.names = F, quote = F, sep = "\t")
