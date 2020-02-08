tmp <- read.delim("~/farm/md5checksum_failed.txt", header = F)
tmp$V1 <- gsub(": FAILED", "", tmp$V1)
tmp$V1 <- gsub(" open or read", "", tmp$V1)
write.table(tmp, "~/farm/failed.txt", quote = F, col.names = F, row.names = F)
