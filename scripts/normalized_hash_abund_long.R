## Get file names
files <- snakemake@input
files <- unlist(files, use.names=FALSE)

## read files into list with each row labelled by sample. 
## normalize files by number of hashes when importing.
ibd_long <- list()
for(i in 1:length(files)){
  print(i)
  sig <- read.csv(files[i])                # read in signature csv as df
  num_hashes <- nrow(sig)                  # get number of hashes
  sig[ , 2] <- sig[ , 2] / num_hashes      # normalize by number of hashes
  sample <- colnames(sig)[2]               # set lib name using sample
  sample <- gsub("_filt.sig", "", sample)  # remove file suffix
  sample <- gsub("^X", "", sample)         # remove X appended by R on import
  sig$sample <- sample                     # set libname as col
  colnames(sig) <- c("minhash", "abund", "sample")
  ibd_long[[i]] <- sig
}

## bind into one dataframe
ibd_long <- do.call(rbind, ibd_long)
write.csv(ibd_long, snakemake@output[['csv']],
          quote = F, row.names = F)

