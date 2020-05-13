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

ibd <- read_feather(snakemake@input[['feather']]) # read in hash abund table
ibd <- as.data.frame(ibd)                         # transform to dataframe
ibd$sample <- gsub("_filt_named\\.sig", "", ibd$sample)
rownames(ibd) <- ibd$sample                       # set sample as rownames
ibd <- ibd[ , -ncol(ibd)]                         # remove the samples column

## read in study metadata
## collapse duplicate libraries so each sample only has one row
info <- read_tsv(snakemake@input[['info']]) %>%
  select(study_accession, library_name, diagnosis) %>%
  mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  filter(library_name %in% rownames(ibd)) %>%
  distinct()

## set validation cohort and remove it from variable selection
# info_validation <- info %>%
#   filter(study_accession == snakemake@params[["validation_study"]]) %>%
#   mutate(library_name = gsub("-", "\\.", library_name))
# ibd_validation <- ibd[rownames(ibd) %in% info_validation$library_name, ]

## remove validation cohort from variable selection
info_novalidation <- info %>%
  filter(study_accession != snakemake@params[["validation_study"]]) %>%
  mutate(library_name = gsub("-", "\\.", library_name))
ibd_novalidation <- ibd[rownames(ibd) %in% info_novalidation$library_name, ]

## make classification vector
## match order of to ibd
info_novalidation <- info_novalidation[match(rownames(ibd_novalidation), info_novalidation$library_name), ]
## make diagnosis var
diagnosis_novalidation <- info_novalidation$diagnosis

# run vita ----------------------------------------------------------------

## perform variant selection
## var.sel.vita calculates p-values based on the empirical null distribution
## from non-positive VIMs as described in Janitza et al. (2015).
ibd_vita <- var.sel.vita(x = ibd_novalidation, y = diagnosis_novalidation, p.t = 0.05,
                        ntree = 10000, mtry.prop = 0.2, nodesize.prop = 0.1,
                        no.threads = snakemake@params[["threads"]], 
                        method = "ranger", type = "classification")
saveRDS(ibd_vita, snakemake@output[["vita_rf"]])

# write files -------------------------------------------------------------

## write predictive hashes
var <- ibd_vita$var                 # separate out selected predictive hashes
var <- gsub("X", "", var)           # remove the X from the beginning of hashes
write.table(var, snakemake@output[['vita_vars']],
            quote = F, col.names = F, row.names = F)

## filter to predictive hashes and write to fie
ibd_filt <- ibd[ , colnames(ibd) %in% var] # subset ibd to hashes in ibd_vita
write.csv(ibd_filt, snakemake@output[['ibd_filt']], quote = F)
