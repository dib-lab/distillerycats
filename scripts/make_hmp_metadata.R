library(dplyr)
library(readr)
hmp <- read.csv("sandbox/08_hmp_validation/hash_abund_tbl/hmp_wide_ibd.csv")
rownames(hmp) <- hmp$sample
hmp <- select(hmp, -sample, -sample.1)
hmp <- as.data.frame(hmp)

info_hmp <- read.csv("~/github/cosmo-kmers/inputs/hmp2_metadata.csv", stringsAsFactors = F) %>%
  filter(data_type == "metagenomics") %>%
  filter(External.ID %in% rownames(hmp))

write_tsv(info_hmp, "inputs/hmp2_mgx_metadata.tsv")
