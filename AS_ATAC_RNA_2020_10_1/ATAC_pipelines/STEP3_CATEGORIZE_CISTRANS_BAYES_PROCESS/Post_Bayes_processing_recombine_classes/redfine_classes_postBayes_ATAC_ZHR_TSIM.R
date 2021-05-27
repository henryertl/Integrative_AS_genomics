#### Goal: Compile everything and run together and CPM transform
library(plyr)
library(dplyr)

rm(list = ls())

setwd("/Users/wittkopp_member/Code")
setwd("/Users/henryertl/Documents/Devs")

##### 1000bp centered intra and inter #####
start <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txStart.txt", header = T)
start$class <- "start"
start$paste_locus <- paste(start[,1], start[,2], start[,3], sep = "_")
start[,c(1:3)] <- NULL

end <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txEnd.txt", header = T)
end$class <- "end"
end$paste_locus <- paste(end[,1], end[,2], end[,3], sep = "_")
end[,c(1:3)] <- NULL

inter <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_TSIM_ATAC_20min1000max_intergenic_NO_TxStart_End.bed", header = F)
inter$paste_locus <- paste(inter$V1, inter$V2, inter$V3, sep = "_")
inter$class <- "inter"
inter[,c(1:3)] <- NULL

intra <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_TSIM_ATAC_20min1000max_intragenic_NO_TxStart_End.bed", header = F)
intra$paste_locus <- paste(intra$V1, intra$V2, intra$V3, sep = "_")
intra$class <- "intra"
intra[,c(1:3)] <- NULL

df_ATAC <- rbind(start, end, inter, intra)

colnames(df_ATAC) <- c("Paste_locus", "class")

# read in full results output
results <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_TSIM_ATAC_20min1000max_no_classes.txt", header = T) %>% unique()
results$Direction <- NULL

results_classes <- left_join(results, df_ATAC, by = "Paste_locus") %>% na.omit() %>% unique()
write.table(results_classes, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_TSIM_ATAC_20min1000max_ALL_classes.txt", sep = "\t", row.names = F, quote = F)
