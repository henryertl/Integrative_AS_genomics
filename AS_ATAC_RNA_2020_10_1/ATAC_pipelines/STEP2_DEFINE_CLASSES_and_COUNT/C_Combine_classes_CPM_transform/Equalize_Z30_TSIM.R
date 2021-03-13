Z30 <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Raw_datatables/ZHR_Z30_ATAC_counts_ALLclasses_centered1000.txt", header = T) %>% unique()
TSIM <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Raw_datatables/ZHR_TSIM_ATAC_counts_ALLclasses_centered1000.txt", header = T) %>% unique()

colSums(Z30[,c(3:ncol(Z30)-1)], na.rm = TRUE)
colSums(TSIM[,c(3:ncol(Z30)-1)], na.rm = TRUE)

# remove columns w/ < 5563321 read sum (mean of Z30)
TSIM <- TSIM[-c(11,12,14,15,16)]

write.table(TSIM, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Raw_datatables/ZHR_TSIM_ATAC_counts_ALLclasses_centered1000_omit_lessZ30MEAN.txt", sep = "\t", quote = F, row.names = F)
