library(readr)
library(dplyr)
library(vegan)

# read in comp and format  ----------------------------------------------
comp <- read_csv(snakemake@input[['comp']]) # read in mat
colnames(comp) <- gsub("_filt", "", colnames(comp))
rownames(comp) <-colnames(comp) 

# read in and formate metadata -------------------------------------------

info <- read_tsv(snakemake@input[['info']]) %>%
  filter(library_name %in% colnames(comp)) %>%
  group_by(library_name, study_accession, diagnosis, subject) %>%
  summarise(read_count = sum(read_count)) %>%
  as.data.frame()

sigs <- read_csv(snakemake@input[['sig_info']]) %>%
  mutate(name = gsub("_filt", "", name)) %>%
  select(name, n_hashes)

info <- left_join(info, sigs, by = c("library_name" = "name"))

# dist and permanova ------------------------------------------------------

dist <- dist(comp)                                       # compute distances
info <- info[match(colnames(comp), info$library_name), ] # sort info by colnames
info$diagnosis <- as.factor(info$diagnosis)              # set factors for model
info$study_accession <- as.factor(info$study_accession)
info$subject <- as.factor(info$subject)

perm <- adonis(dist ~ diagnosis + study_accession + read_count + n_hashes, 
               data = info, 
               permutations = 10000)

write.csv(as.data.frame(perm$aov.tab), snakemake@output[['perm']], quote = F)


