library(ranger)
library(dplyr)
library(readr)
library(feather)
library(Pomona)

set.seed(1)

# perform variable selection within random forests on minhashes (k 31,
# scaled 2000, with hashes that were present only once across all samples
# removed). Uses the Pomona package vita implementation, which wraps the ranger
# package. Saves output as RDS for faster loading of output objects into
# subsequent R sessions.

## read in data ------------------------------------------------------------

## format hash table (samples x features)

hashes <- read_feather(snakemake@input[['feather']]) # read in hash abund table
hashes <- as.data.frame(hashes)                      # transform to dataframe
hashes$sample <- gsub("_filt_named\\.sig", "", hashes$sample)
rownames(hashes) <- hashes$sample                    # set sample as rownames
hashes <- hashes[ , -ncol(hashes)]                   # remove the samples column

## read in study metadata
## change "-" to "." in sample names, as some programs do this automatically
## collapse duplicate samples so each sample only has one row
info <- read_tsv(snakemake@input[['info']]) %>%
  mutate(sample = gsub("\\-", "\\.", sample)) %>%
  filter(sample %in% rownames(hashes)) %>%
  distinct()

## remove validation cohort from variable selection
info_novalidation <- info %>%
  filter(study_accession != snakemake@params[["validation_study"]]) %>%
  mutate(sample = gsub("-", "\\.", sample))
hashes_novalidation <- hashes[rownames(hashes) %in% info_novalidation$sample, ]

## make classification vector
## match order of to hashes
info_novalidation <- info_novalidation[match(rownames(hashes_novalidation), info_novalidation$sample), ]
## make var var
var_novalidation <- info_novalidation$var

# run vita ----------------------------------------------------------------

## perform variant selection
## var.sel.vita calculates p-values based on the empirical null distribution
## from non-positive VIMs as described in Janitza et al. (2015).
var_vita <- var.sel.vita(x = hashes_novalidation, y = var_novalidation, p.t = 0.05,
                         ntree = 10000, mtry.prop = 0.2, nodesize.prop = 0.1,
                         no.threads = snakemake@params[["threads"]], 
                         method = "ranger", type = "classification")
saveRDS(var_vita, snakemake@output[["vita_rf"]])

# write files -------------------------------------------------------------

## write predictive hashes
var <- var_vita$var                 # separate out selected predictive hashes
var <- gsub("X", "", var)           # remove the X from the beginning of hashes
write.table(var, snakemake@output[['vita_vars']],
            quote = F, col.names = F, row.names = F)

## filter to predictive hashes and write to fie
hashes_filt <- hashes[ , colnames(hashes) %in% var] # subset ibd to hashes in ibd_vita
write.csv(hashes_filt, snakemake@output[['var_filt']], quote = F)
