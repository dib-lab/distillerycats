library(readr)
library(dplyr)
library(ggplot2)
library(ranger)
library(mlr)
library(tuneRanger)
source(snakemake@input[['eval_model']])
source(snakemake@input[['ggconfusion']])

set.seed(1)

var_filt <- read_csv(snakemake@input[['var_filt']])
var_filt <- as.data.frame(var_filt)
rownames(var_filt) <- var_filt$X1
var_filt <- var_filt[ , -1]

## read in study metadata
## change "-" in sample names to "." as some programs automatically make this change
## filter to contain only samples that are in var_filt
## only keep distinct rows.
info <- read_tsv(snakemake@input[['info']]) %>%
  mutate(sample = gsub("\\-", "\\.", sample)) %>%
  filter(sample %in% rownames(var_filt)) %>%
  distinct()

## set validation cohort and remove it from variable selection
info_validation <- info %>%
  filter(study_accession == snakemake@params[['validation_study']]) %>%
  mutate(sample = gsub("-", "\\.", sample))
var_validation <- var_filt[rownames(var_filt) %in% info_validation$sample, ]
# match order of to var_filt
info_validation <- info_validation[order(match(info_validation$sample, rownames(var_validation))), ]
# check names
all.equal(info_validation$sample, rownames(var_validation))
# make var col
var_validation <- info_validation$var


## remove validation cohort from training data
# using tuneRanger, we do not need to use a train/test/validation framework.
# Instead, tuneRanger does not need a test set because each tree is only trained 
# on a subset of the data (bag), so we can use the rest (out of bag) to obtain 
# an unbiased performance estimation of a single tree and therefore of all trees.
# see: https://github.com/PhilippPro/tuneRanger/issues/8

info_novalidation <- info %>%
  filter(study_accession != snakemake@params[['validation_study']]) %>%
  mutate(sample = gsub("-", "\\.", sample))
var_novalidation <- var_filt[rownames(var_filt) %in% info_novalidation$sample, ]
# match order of to var_filt
info_novalidation <- info_novalidation[order(match(info_novalidation$sample, rownames(var_novalidation))), ]
# check names
all.equal(info_novalidation$sample, rownames(var_novalidation))
# make var col
var_novalidation <- info_novalidation$var

# Include classification vars as cols in df
var_novalidation$var <- var_novalidation
var_validation$var <- var_validation

# tune ranger -------------------------------------------------------------

# Make an mlr task with the var_train dataset here 
tmp <- var_novalidation
colnames(tmp) <-  make.names(colnames(tmp))
var_task <- makeClassifTask(data = tmp, target = "var")
# Run tuning process
res <- tuneRanger(var_task, num.threads = snakemake@params[['threads']])

# write model parameters to a file
write_tsv(res$recommended.pars, snakemake@output[['recommended_pars']])

# build optimal model ----------------------------------------------------------

# extract model parameters and use to build an optimal RF

# use model parameters to build optimized RF
var_novalidation$var <- as.factor(var_novalidation$var)
optimal_rf <- ranger(
  dependent.variable.name = "var",
  mtry            = res$recommended.pars$mtry,
  num.trees       = 10000,
  data            = var_novalidation,
  sample.fraction = res$recommended.pars$sample.fraction,
  min.node.size   = res$recommended.pars$min.node.size,
  seed            = 1,
  importance      = 'impurity'
)

saveRDS(optimal_rf, file = snakemake@output[['optimal_rf']])

# evaluate the accuracy of the model and generate a confusion matrix
# training data
evaluate_model(optimal_ranger = optimal_rf, 
               data = var_novalidation, 
               reference_class = var_novalidation, 
               set = "novalidation", 
               study_as_validation = snakemake@params[['validation_study']],
               accuracy_csv = snakemake@output[['training_accuracy']],
               confusion_pdf = snakemake@output[['training_confusion']])
# validation data
evaluate_model(optimal_ranger = optimal_rf, 
               data = var_validation, 
               reference_class = var_validation, 
               set = "validation", 
               study_as_validation = snakemake@params[['validation_study']],
               accuracy_csv = snakemake@output[['validation_accuracy']],
               confusion_pdf = snakemake@output[['validation_confusion']])
