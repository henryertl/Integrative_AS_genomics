#### Goal: Compile everything and run together and CPM transform
library(plyr)
library(dplyr)

rm(list = ls())

setwd("/Users/wittkopp_member/Code")

##### 1000bp centered intra and inter #####
start <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_txStart500_counts.bed", header = F)
start$V4 <- NULL
start$class <- "start"
end <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_txEnd500_counts.bed", header = F)
end$V4 <- NULL
end$class <- "end"
inter <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_intergenic_peak_counts_center1000.bed", header = F)
inter$class <- "inter"
intra <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_ATAC_intragenic_peak_counts_center1000.bed", header = F)
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

df_ATAC_class_locus_key <- df_ATAC[,c("class", "Paste_locus")]

# read in full results output
results <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq/Data_tables/Full_results_output_ZHR_Z30_ATAC_20min_centered.txt", header = T)
results$Paste_locus <- paste(results$chrom, results$start, results$end, sep = "_")
results$Direction <- NULL

test <- left_join(results, df_ATAC_class_locus_key, by = "Paste_locus") %>% na.omit() %>% unique()
