library(readr)
library(dplyr)
library(ggplot2)
library(ranger)
library(mlr)
library(tuneRanger)
source(snakemake@input[['eval_model']])
source(snakemake@input[['ggconfusion']])

set.seed(1)

ibd_filt <- read_csv(snakemake@input[['ibd_filt']])
ibd_filt <- as.data.frame(ibd_filt)
rownames(ibd_filt) <- ibd_filt$X1
ibd_filt <- ibd_filt[ , -1]

## read in study metadata
## collapse duplicate libraries so each sample only has one row
info <- read_tsv(snakemake@input[['info']]) %>%
  select(study_accession, library_name, diagnosis) %>%
  mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  filter(library_name %in% rownames(ibd_filt)) %>%
  distinct()

## set validation cohort and remove it from variable selection
info_validation <- info %>%
  filter(study_accession == snakemake@params[['validation_study']]) %>%
  mutate(library_name = gsub("-", "\\.", library_name))
ibd_validation <- ibd_filt[rownames(ibd_filt) %in% info_validation$library_name, ]
# match order of to ibd_filt
info_validation <- info_validation[order(match(info_validation$library_name, rownames(ibd_validation))), ]
# check names
all.equal(info_validation$library_name, rownames(ibd_validation))
# make diagnosis var
diagnosis_validation <- info_validation$diagnosis


## remove validation cohort from training data
# using tuneRanger, we do not need to use a train/test/validation framework.
# Instead, tuneRanger does not need a test set because each tree is only trained 
# on a subset of the data (bag), so we can use the rest (out of bag) to obtain 
# an unbiased performance estimation of a single tree and therefore of all trees.
# see: https://github.com/PhilippPro/tuneRanger/issues/8

info_novalidation <- info %>%
  filter(study_accession != snakemake@params[['validation_study']]) %>%
  mutate(library_name = gsub("-", "\\.", library_name))
ibd_novalidation <- ibd_filt[rownames(ibd_filt) %in% info_novalidation$library_name, ]
# match order of to ibd_filt
info_novalidation <- info_novalidation[order(match(info_novalidation$library_name, rownames(ibd_novalidation))), ]
# check names
all.equal(info_novalidation$library_name, rownames(ibd_novalidation))
# make diagnosis var
diagnosis_novalidation <- info_novalidation$diagnosis

# Include classification vars as cols in df
ibd_novalidation$diagnosis <- diagnosis_novalidation
ibd_validation$diagnosis <- diagnosis_validation

# tune ranger -------------------------------------------------------------

# Make an mlr task with the ibd_train dataset here 
tmp <- ibd_novalidation
colnames(tmp) <-  make.names(colnames(tmp))
ibd_task <- makeClassifTask(data = tmp, target = "diagnosis")
# Run tuning process
res <- tuneRanger(ibd_task, num.threads = snakemake@params[['threads']])

# write model parameters to a file
write_tsv(res$recommended.pars, snakemake@output[['recommended_pars']])

# build optimal model ----------------------------------------------------------

# extract model parameters and use to build an optimal RF

# use model parameters to build optimized RF
ibd_novalidation$diagnosis <- as.factor(ibd_novalidation$diagnosis)
optimal_rf <- ranger(
  dependent.variable.name = "diagnosis",
  mtry            = res$recommended.pars$mtry,
  num.trees       = 10000,
  data            = ibd_novalidation,
  sample.fraction = res$recommended.pars$sample.fraction,
  min.node.size   = res$recommended.pars$min.node.size,
  seed            = 1,
  importance      = 'impurity'
)

saveRDS(optimal_rf, file = snakemake@output[['optimal_rf']])

# evaluate the accuracy of the model and generate a confusion matrix
# training data
evaluate_model(optimal_ranger = optimal_rf, 
               data = ibd_novalidation, 
               reference_class = diagnosis_novalidation, 
               set = "novalidation", 
               study_as_validation = snakemake@params[['validation_study']],
               accuracy_csv = snakemake@output[['training_accuracy']],
               confusion_pdf = snakemake@output[['training_confusion']])
# validation data
evaluate_model(optimal_ranger = optimal_rf, 
               data = ibd_validation, 
               reference_class = diagnosis_validation, 
               set = "validation", 
               study_as_validation = snakemake@params[['validation_study']],
               accuracy_csv = snakemake@output[['validation_accuracy']],
               confusion_pdf = snakemake@output[['validation_confusion']])
