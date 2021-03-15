#### Goal: Compile everything and run together and CPM transform
library(plyr)
library(dplyr)

rm(list = ls())

setwd("/Users/wittkopp_member/Code")

##### 1000bp centered intra and inter #####
start <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_TSIM_ATAC_txStart500_counts.bed", header = F) %>% unique()
start$V4 <- NULL
start$class <- "start"
end <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_TSIM_ATAC_txEnd500_counts.bed", header = F) %>% unique()
end$V4 <- NULL
end$class <- "end"
inter <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_TSIM_ATAC_intergenic_peak_counts_center1000.bed", header = F) %>% unique()
inter$class <- "inter"
intra <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_TSIM_ATAC_intragenic_peak_counts_center1000.bed", header = F) %>% unique()
intra$class <- "intra"

colnames(start) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "class")
colnames(end) <- c("chrom", "start", "end",  "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "class")
colnames(inter) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "class")
colnames(intra) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "class")

df_ATAC <- rbind(start, end, inter, intra)
df_ATAC$Paste_locus <- paste(df_ATAC$chrom, df_ATAC$start, df_ATAC$end, sep = "_")

# generate paste_locus key to use to assign downsamp data
classes_key_TSIM <- df_ATAC[,c(1:3, 16)]
classes_key_TSIM$Paste_locus <- paste(classes_key_TSIM$chrom, classes_key_TSIM$start, classes_key_TSIM$end, sep = "_")
classes_key_TSIM[,c(1:3)] <- NULL

## reformat downsampled dataset to have classes and then append to rest of data
down_samp <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Raw_datatables/ZHR_TSIM_combined_ALLclasses_ATAC_downsampled_d_merge8_Z30mean_sub.txt", header = T)
down_samp$Paste_locus <- paste(down_samp$chrom, down_samp$start, down_samp$end, sep = "_")
down_samp_classes <- left_join(down_samp, classes_key_TSIM, by = "Paste_locus")

df_ATAC_overlap_down_samp <- join_all(list(down_samp_classes, df_ATAC), type = "full", by = "Paste_locus") %>% na.omit() %>% unique()

df_ATAC_final <- df_ATAC_overlap_down_samp[,c(10:12,9:3,23,24,2,25,26,22)]

colnames(df_ATAC_final) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "class")

write.table(df_ATAC_final, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Raw_datatables/ZHR_TSIM_ATAC_counts_ALLclasses_centered1000_downsamp_nondownsamp_comb.txt", row.names = F, quote = F, sep = "\t")

full_dataset <- df_ATAC_final %>% as.data.frame() %>% na.omit() %>% unique()

full_dataset <- full_dataset[full_dataset$P1_1 >= 20 & full_dataset$P1_2 >= 20 & full_dataset$P1_3 >= 20 & full_dataset$P2_1 >= 20 & full_dataset$P2_2 >= 20 &
full_dataset$P2_3 >= 20 & full_dataset$HYB_1_P1 + full_dataset$HYB_1_P2 >= 20 & full_dataset$HYB_2_P1 + full_dataset$HYB_2_P2 >= 20 & full_dataset$HYB_3_P1 + full_dataset$HYB_3_P2 >= 20,]

df_ATAC <- full_dataset %>% as.data.frame()

write.table(df_ATAC, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Raw_datatables/ZHR_TSIM_ATAC_counts_ALLclasses_downsamp_nondownsamp_comb_20min_centered1000.txt", row.names = F, quote = F, sep = "\t")

# CPM transformation
df_ATAC_CPM <- matrix(ncol = ncol(df_ATAC), nrow = nrow(df_ATAC)) %>% as.data.frame()
df_ATAC_CPM[,ncol(df_ATAC)] <- df_ATAC$class

for(i in 4:ncol(df_ATAC)-1) {
  df_ATAC_CPM[,i] <- (df_ATAC[,i]/sum(df_ATAC[,i]))*1000000
}

df_ATAC_CPM$V1 <- df_ATAC$chrom
df_ATAC_CPM$V2 <- df_ATAC$start
df_ATAC_CPM$V3 <- df_ATAC$end
colnames(df_ATAC_CPM) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1","HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "class")

## output == CPM transformed values for each sample with class assignments
### 11626 txStart
### 11619 txEnd
### 2098 inter
### 4801 intra

write.table(df_ATAC_CPM, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/CPM_transformed_datatables/ZHR_TSIM_ATAC_counts_downsamp_nondownsamp_comb_ALLclasses_20min_CPM_centered1000.txt", row.names = F, quote = F, sep = "\t")

# also create dataframe with only regions also analyzed in ZHR - Z30
Z30 <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/CPM_transformed_datatables/ZHR_Z30_ATAC_counts_ALLclasses_20min_CPM_centered1000.txt", header = T)
colnames(Z30) <- c("chromZ30", "startZ30", "endZ30", "P1_1Z30", "P1_2Z30", "P1_3Z30", "P2_1Z30", "P2_2Z30", "P2_3Z30",
"HYB_1_P1Z30", "HYB_2_P1Z30", "HYB_3_P1Z30","HYB_1_P2Z30", "HYB_2_P2Z30", "HYB_3_P2Z30", "classZ30")
Z30$Paste_locus <- paste(Z30$chrom, Z30$start, Z30$end, sep = "_")

df_ATAC_CPM$Paste_locus <- paste(df_ATAC_CPM$chrom, df_ATAC_CPM$start, df_ATAC_CPM$end, sep = "_")

TSIM_only_Z30overlap <- left_join(df_ATAC_CPM, Z30, by = "Paste_locus", type = "full") %>% na.omit() %>% unique()
TSIM_only_Z30overlap <- TSIM_only_Z30overlap[,c(1:16)]

write.table(TSIM_only_Z30overlap, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/CPM_transformed_datatables/ZHR_TSIM_ATAC_counts_downsamp_nondownsamp_comb_ALLclasses_20min_CPM_centered1000_onlyZ30overlap.txt", row.names = F, quote = F, sep = "\t")
