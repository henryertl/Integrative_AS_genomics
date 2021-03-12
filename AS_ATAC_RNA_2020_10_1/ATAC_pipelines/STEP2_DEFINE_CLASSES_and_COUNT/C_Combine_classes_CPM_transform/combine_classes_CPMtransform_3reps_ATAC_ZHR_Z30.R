#### Goal: Compile everything and run together and CPM transform
library(plyr)
library(dplyr)

rm(list = ls())

setwd("/Users/wittkopp_member/Code")

##### 1000bp centered intra and inter #####
start <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_txStart500_counts.bed", header = F)
start$V4 <- NULL
end <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_txEnd500_counts.bed", header = F)
end$V4 <- NULL
inter <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_intergenic_peak_counts_center1000.bed", header = F)
intra <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_intragenic_peak_counts_center1000.bed", header = F)

colnames(start) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")
colnames(end) <- c("chrom", "start", "end",  "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")
colnames(inter) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")
colnames(intra) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")

df_ATAC <- rbind(start, end, inter, intra)

full_dataset <- df_ATAC %>% as.data.frame() %>% na.omit() %>% unique()

full_dataset <- full_dataset[full_dataset$P1_1 >= 20 & full_dataset$P1_2 >= 20 & full_dataset$P1_3 >= 20 & full_dataset$P2_1 >= 20 & full_dataset$P2_2 >= 20 &
full_dataset$P2_3 >= 20 & full_dataset$HYB_1_P1 + full_dataset$HYB_1_P2 >= 20 & full_dataset$HYB_2_P1 + full_dataset$HYB_2_P2 >= 20 & full_dataset$HYB_3_P1 + full_dataset$HYB_3_P2 >= 20,]

df_ATAC <- full_dataset %>% as.data.frame()

# CPM transformation
df_ATAC_CPM <- matrix(ncol = ncol(df_ATAC) - 3, nrow = nrow(df_ATAC)) %>% as.data.frame()

for(i in 4:ncol(df_ATAC)) {
  df_ATAC_CPM[,i] <- (df_ATAC[,i]/sum(df_ATAC[,i]))*1000000
}

df_ATAC_CPM$V1 <- df_ATAC$chrom
df_ATAC_CPM$V2 <- df_ATAC$start
df_ATAC_CPM$V3 <- df_ATAC$end
colnames(df_ATAC_CPM) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1","HYB_1_P2", "HYB_2_P2", "HYB_3_P2")


write.table(df_ATAC_CPM, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq/Data_tables/ZHR_Z30_ATAC_counts_ALLclasses_20min_CPM_centered1000.txt", row.names = F, quote = F, sep = "\t")


##### normal intra and inter #####
start <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_txStart500_counts.bed", header = F)
start$V4 <- NULL
end <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_txStart500_counts.bed", header = F)
end$V4 <- NULL
inter <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_intergenic_peak_counts.bed", header = F)
intra <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_intragenic_peak_counts.bed", header = F)

colnames(start) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")
colnames(end) <- c("chrom", "start", "end",  "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")
colnames(inter) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")
colnames(intra) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")

df_ATAC <- rbind(start, end, inter, intra)

full_dataset <- df_ATAC %>% as.data.frame() %>% na.omit()

full_dataset <- full_dataset[full_dataset$P1_1 >= 20 & full_dataset$P1_2 >= 20 & full_dataset$P1_3 >= 20 & full_dataset$P2_1 >= 20 & full_dataset$P2_2 >= 20 &
full_dataset$P2_3 >= 20 & full_dataset$HYB_1_P1 + full_dataset$HYB_1_P2 >= 20 & full_dataset$HYB_2_P1 + full_dataset$HYB_2_P2 >= 20 & full_dataset$HYB_3_P1 + full_dataset$HYB_3_P2 >= 20,]

df_ATAC <- full_dataset %>% as.data.frame()

# CPM transformation
df_ATAC_CPM <- matrix(ncol = ncol(df_ATAC) - 3, nrow = nrow(df_ATAC)) %>% as.data.frame()

for(i in 4:ncol(df_ATAC)) {
  df_ATAC_CPM[,i] <- (df_ATAC[,i]/sum(df_ATAC[,i]))*1000000
}

df_ATAC_CPM$V1 <- df_ATAC$chrom
df_ATAC_CPM$V2 <- df_ATAC$start
df_ATAC_CPM$V3 <- df_ATAC$end
colnames(df_ATAC_CPM) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1","HYB_1_P2", "HYB_2_P2", "HYB_3_P2")
write.table(df_ATAC_CPM, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq/Data_tables/ZHR_Z30_ATAC_counts_ALLclasses_20min_CPM.txt", row.names = F, quote = F, sep = "\t")
