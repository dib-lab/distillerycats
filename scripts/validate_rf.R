
library(dplyr)
library(readr)
library(ranger)

set.seed(1)

# validation by study -----------------------------------------------------

optimal_ranger <- readRDS(snakemake@input[['optimal_rf']])

# or read in subsetted ibd validation filter file
ibd_validation_filt <- read.csv(snakemake@input[['ibd_validation']],
                                row.names = 1)

info_validation <- read_tsv(snakemake@input[['info']]) %>% 
  select(study_accession, library_name, diagnosis) %>%
  mutate(library_name = gsub("\\-", "\\.", library_name)) %>%
  filter(library_name %in% rownames(ibd_validation_filt)) %>%
  distinct()


# SRP057027 ----------------------------------------------------------------

info_srp <- filter(info_validation, study_accession == "SRP057027")
ibd_srp <- ibd_validation_filt[rownames(ibd_validation_filt) %in% info_srp$library_name, ]
info_srp <- info_srp[match(rownames(ibd_srp), info_srp$library_name), ]
ibd_srp$diagnosis <- info_srp$diagnosis

pred_srp <- predict(optimal_ranger, ibd_srp)
pred_srp_tab <- table(observed = ibd_srp$diagnosis, predicted = pred_srp$predictions)
write.table(pred_srp_tab, snakemake@output[['pred_srp']])

pred_srp_df <- data.frame(sample = rownames(ibd_srp),
                            diagnosis = ibd_srp$diagnosis,
                            prediction = pred_srp$predictions)

write.csv(pred_srp_df, snakemake@output[['pred_srp_df']],
          quote = F, row.names = F)
# PRJNA285949 -------------------------------------------------------------

info_prjna <- filter(info_validation, study_accession == "PRJNA385949")
ibd_prjna <- ibd_validation_filt[rownames(ibd_validation_filt) %in% info_prjna$library_name, ]
info_prjna <- info_prjna[match(rownames(ibd_prjna), info_prjna$library_name), ]
ibd_prjna$diagnosis <- info_prjna$diagnosis

pred_prjna <- predict(optimal_ranger, ibd_prjna)
pred_prjna_tab <- table(observed = ibd_prjna$diagnosis, predicted = pred_prjna$predictions)
write.table(pred_srp_tab, snakemake@output[['pred_prjna']])

pred_prjna_df <- data.frame(sample = rownames(ibd_prjna),
                            diagnosis = ibd_prjna$diagnosis,
                            prediction = pred_prjna$predictions)
write.csv(pred_prjna_df, snakemake@output[['pred_prjna_df']],
          quote = F, row.names = F)


