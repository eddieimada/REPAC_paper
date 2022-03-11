library(tidyverse)

# This scripts creates the base expression profiles for the 2 conditions simulated with simFlux. The base profiles are provided in this repository.
## read annotation
anno.sel <- read_rds("anno_sel.rds")
gns <- unique(anno.sel$gene_id)
# Load expression profiles for the simulation
proA <- read_tsv("../data/parA.pro", col_names = F)
set.seed(777)
# Randomly assign a ratio of long:short isoform between 10:90 to 35:65

idx <- sample(c(1,2), length(gns), replace = T)
truth <- cbind(gns,idx)
pct1 <- unlist(map(idx, ~{
    if (.x == 1) {
        r <- sample(seq(0.1,0.35, by=0.01), 1)
        c(1-r,r)
    } else {
        r <- sample(seq(0.1,0.35, by=0.01), 1)
        c(r,1-r)
    }
}))

proA <- bind_cols(proA,pct=pct1)
proA <- proA %>%
    mutate(id=gsub("L|S", "", X2)) %>%
    group_by(id) %>%
    mutate(X6=round(sum(X6)*pct)) %>%
    ungroup() %>%
    mutate(X5=X6/sum(X6)) %>%
    select(-pct,-id)
    
# Same for the second condition
proB <- read_tsv("/Users/eimada/Downloads/flux-simulator-1.2.1/bin/parB.pro", col_names = F)

pct2 <- unlist(map(idx, ~{
    if (.x == 2) {
        r <- sample(seq(0.1,0.35, by=0.01), 1)
        c(1-r,r)
    } else {
        r <- sample(seq(0.1,0.35, by=0.01), 1)
        c(r,1-r)
    }
}))

proB <- bind_cols(proB,pct=pct2)
proB <- proB %>%
    mutate(id=gsub("L|S", "", X2)) %>%
    group_by(id) %>%
    mutate(X6=round(sum(X6)*pct)) %>%
    ungroup() %>%
    mutate(X5=X6/sum(X6)) %>%
    select(-pct,-id)

# Save profiles to use with simFlux
write.table(proA, file = "simA_init2iso.pro",col.names = F, row.names = F, quote = F, sep = "\t")
write.table(proB, file = "simB_init2iso.pro",col.names = F, row.names = F, quote = F, sep = "\t")
