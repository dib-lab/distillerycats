library(readr)
library(dplyr)
library(vegan)

# read in comp and format  ----------------------------------------------
comp <- read_csv(snakemake@input[['comp']]) # read in mat
colnames(comp) <- gsub("_filt", "", colnames(comp))
rownames(comp) <-colnames(comp) 

# read in and formate metadata -------------------------------------------

info <- read_tsv(snakemake@input[['info']]) %>%
  filter(sample %in% colnames(comp)) %>%
  as.data.frame()

sigs <- read_csv(snakemake@input[['sig_info']]) %>%
  mutate(name = gsub("_filt", "", name)) %>%
  select(name, n_hashes)

info <- left_join(info, sigs, by = c("sample" = "name"))

# dist and permanova ------------------------------------------------------

dist <- dist(comp)                                 # compute distances
info <- info[match(colnames(comp), info$sample), ] # sort info by colnames
info$var <- as.factor(info$var)                    # set factors for model
info$study <- as.factor(info$study)


perm <- adonis(dist ~ var + study + n_hashes, 
               data = info, 
               permutations = 10000)

write.csv(as.data.frame(perm$aov.tab), snakemake@output[['perm']], quote = F)