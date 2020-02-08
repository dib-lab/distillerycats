library(ranger)
library(dplyr)
library(readr)
library(feather)
library(Pomona)

set.seed(1)

# perform variable selection within random forests on ibd minhashes (k 31, 
# scaled 2000, with hashes that were present only once across all samples
# removed). Uses the Pomona package vita implementation, which wraps the ranger 
# package. Saves output as RDS for faster loading of output objects into 
# subsequent R sessions. 


## read in data ------------------------------------------------------------

## format hash table (samples x features)

ibd <- read_feather(snakemake@input[['feather']]) # read in hash abund table
ibd <- as.data.frame(ibd)                         # transform to dataframe
rownames(ibd) <- ibd$sample                       # set sample as rownames
ibd <- ibd[ , -ncol(ibd)]                         # remove the samples column

## read in study metadata
info <- read_tsv(snakemake@input[['info']]) %>% 
  select(study_accession, library_name, diagnosis) %>%
  mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  filter(library_name %in% rownames(ibd)) %>%
  distinct()

## set validation cohorts (e.g. rm time series from training and testing)
info_validation <- info %>% 
  filter(study_accession %in% c("SRP057027", "PRJNA385949"))

info_novalidation <- info %>% 
  filter(!study_accession %in% c("SRP057027", "PRJNA385949"))

## remove validation from ibd
ibd_novalidation <- ibd[rownames(ibd) %in% info_novalidation$library_name, ]

## make classification vector
## match order of classification to ibd
info_novalidation <- info_novalidation[match(rownames(ibd_novalidation), info_novalidation$library_name), ]
## make diagnosis var
diagnosis <- info_novalidation$diagnosis 

## split to test and train
train <- sample(nrow(ibd_novalidation), 0.7*nrow(ibd_novalidation), replace = FALSE)
train_set <- ibd_novalidation[train, ]
test_set <- ibd_novalidation[-train, ]

diagnosis_train <- diagnosis[train]
diagnosis_test <- diagnosis[-train]

# run vita ----------------------------------------------------------------

## perform variant selection
## var.sel.vita calculates p-values based on the empirical null distribution 
## from non-positive VIMs as described in Janitza et al. (2015). 
ibd_vita <- var.sel.vita(x = train_set, y = diagnosis_train, p.t = 0.05, 
                         ntree = 5000, mtry.prop = 0.2, nodesize.prop = 0.1, 
                         no.threads = 10, method = "ranger", 
                         type = "classification")
saveRDS(ibd_vita, snakemake@output[["vita_rf"]])

# write files -------------------------------------------------------------

## write predictive hashes
var <- ibd_vita$var                 # separate out selected predictive hashes
var <- gsub("X", "", var)           # remove the X from the beginning of hashes
write.table(var, snakemake@output[['vita_vars']], 
            quote = F, col.names = F, row.names = F)

## filter to predictive hashes and write training/testing set (novalidation) to files
ibd_filt <- ibd_novalidation[ , colnames(ibd_novalidation) %in% var] # subset ibd to hashes in ibd_vita
write.csv(ibd_filt, snakemake@output[['ibd_novalidation']], quote = F)
write.table(diagnosis, snakemake@output[['ibd_novalidation_diagnosis']], 
            row.names = F, quote = F)

## filter to predictive hashes and write validation set to files
ibd_validation <- ibd[rownames(ibd) %in% info_validation$library_name, ]
ibd_validation_filt <- ibd_validation[ , colnames(ibd_validation) %in% var] # subset ibd to hashes in ibd_vita
write.csv(ibd_validation_filt, snakemake@output[['ibd_validation']], quote = F)

