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

write.table(df_ATAC, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Raw_datatables/ZHR_TSIM_ATAC_counts_ALLclasses_centered1000.txt", row.names = F, quote = F, sep = "\t")

full_dataset <- df_ATAC %>% as.data.frame() %>% na.omit() %>% unique()

full_dataset <- full_dataset[full_dataset$P1_1 >= 20 & full_dataset$P1_2 >= 20 & full_dataset$P1_3 >= 20 & full_dataset$P2_1 >= 20 & full_dataset$P2_2 >= 20 &
full_dataset$P2_3 >= 20 & full_dataset$HYB_1_P1 + full_dataset$HYB_1_P2 >= 20 & full_dataset$HYB_2_P1 + full_dataset$HYB_2_P2 >= 20 & full_dataset$HYB_3_P1 + full_dataset$HYB_3_P2 >= 20,]

df_ATAC <- full_dataset %>% as.data.frame()

write.table(df_ATAC, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Raw_datatables/ZHR_TSIM_ATAC_counts_ALLclasses_20min_centered1000.txt", row.names = F, quote = F, sep = "\t")

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
### 6154 txStart
### 6091 txEnd
### 842 inter
### 1537 intra

write.table(df_ATAC_CPM, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/CPM_transformed_datatables/ZHR_TSIM_ATAC_counts_ALLclasses_20min_CPM_centered1000.txt", row.names = F, quote = F, sep = "\t")
