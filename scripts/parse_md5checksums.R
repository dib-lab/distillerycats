setwd("~/github/ibd")

library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)

files <- list.files("inputs/ENA/", full.names = T)
prj <- files %>%
  map(read_tsv) %>%
  reduce(rbind) %>%
  filter(library_strategy != "AMPLICON") %>%
  filter(library_strategy != "RNA-Seq") %>%
  filter(read_count >= 1000000) %>%
  filter(library_layout  == "PAIRED")

table(prj$study_accession)
# sample_alias has disease status for PRJEB2054
# PRJNA385949 time series (use as validation)
# pediatric crohn's data set: https://raw.githubusercontent.com/louiejtaylor/sbx_lewis2015/master/metadata/metadata_with_SRR.csv
# PRJNA400072 matches library name

md5 <- select(prj, fastq_ftp,  fastq_md5)


# make fastqs as separate --------------------------------------------------

md5 <- separate(data = md5, 
                col = fastq_ftp,
                into = c("fastq_ftp_1", "fastq_ftp_2", "fastq_ftp_3"), 
                sep = ";", remove = F, extra = "merge", fill = "right")
# some fastq files have 3 download links, where the first one contains what is
# probably a merged or interleaved file. Remove these from fastq download links. 
# fastq 1 download:

fastq_1 <- vector()
for(i in 1:length(md5$fastq_ftp)){
  if(grepl(pattern = "_1", x = md5$fastq_ftp_1[i])){
    fastq_1[i] <- md5$fastq_ftp_1[i] 
  } else {
    fastq_1[i] <- md5$fastq_ftp_2[i]
  }
}

fastq_2 <- vector()
for(i in 1:length(md5$fastq_ftp)){
  if(grepl(pattern = "_1", x = md5$fastq_ftp_1[i])){
    fastq_2[i] <- md5$fastq_ftp_2[i] 
  } else {
    fastq_2[i] <- md5$fastq_ftp_3[i]
  }
}

md5$fastq_final_ftp1 <- fastq_1
md5$fastq_final_ftp2 <- fastq_2
md5 <- select(md5, fastq_ftp_1, fastq_final_ftp1, fastq_final_ftp2, fastq_md5)


# separate md5sums too ----------------------------------------------------


md5 <- separate(data = md5, 
                col = fastq_md5,
                into = c("md5_1", "md5_2", "md5_3"), 
                sep = ";", remove = F, extra = "merge", fill = "right")

md5_1 <- vector()
for(i in 1:length(md5$md5_1)){
  if(grepl(pattern = "_1", x = md5$fastq_ftp_1[i])){
    md5_1[i] <- md5$md5_1[i] 
  } else {
    md5_1[i] <- md5$md5_2[i]
  }
}

md5_2 <- vector()
for(i in 1:length(md5$md5_1)){
  if(grepl(pattern = "_1", x = md5$fastq_ftp_1[i])){
    md5_2[i] <- md5$md5_2[i] 
  } else {
    md5_2[i] <- md5$md5_3[i]
  }
}

md5$md5_1 <- md5_1
md5$md5_2 <- md5_2


# make forward and reverse ------------------------------------------------

md5_forward <- select(md5, md5_1, fastq_final_ftp1)
md5_forward$fastq_final_ftp1 <- gsub("ftp.sra.ebi.ac.uk\\/vol1\\/fastq\\/[SERR]{2,3}[0-9]{2,3}\\/[SERR]{2,3}[0-9]{3,7}\\/", 
                                "", md5_forward$fastq_final_ftp1)
md5_forward$fastq_final_ftp1 <- gsub("ftp.sra.ebi.ac.uk\\/vol1\\/fastq\\/[SERR]{2,3}[0-9]{2,3}\\/[0-9]{2,3}\\/[SERR]{2,3}[0-9]{3,7}\\/",
                                "", md5_forward$fastq_final_ftp1)
colnames(md5_forward) <- c('md5sum', "file")


md5_reverse <- select(md5, md5_2, fastq_final_ftp2)
md5_reverse$fastq_final_ftp2 <- gsub("ftp.sra.ebi.ac.uk\\/vol1\\/fastq\\/[SERR]{2,3}[0-9]{2,3}\\/[SERR]{2,3}[0-9]{3,7}\\/", 
                                "", md5_reverse$fastq_final_ftp2)
md5_reverse$fastq_final_ftp2<- gsub("ftp.sra.ebi.ac.uk\\/vol1\\/fastq\\/[SERR]{2,3}[0-9]{2,3}\\/[0-9]{2,3}\\/[SERR]{2,3}[0-9]{3,7}\\/",
                                "", md5_reverse$fastq_final_ftp2)

colnames(md5_reverse) <- c('md5sum', "file")

md5 <- rbind(md5_forward, md5_reverse)

write.table(x = md5, file = "inputs/md5checksums.txt", row.names = F, 
            quote = F, col.names = F, sep = "  ")

