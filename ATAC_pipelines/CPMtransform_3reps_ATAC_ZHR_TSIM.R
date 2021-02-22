#### Goal: Condense exon count to a total genic count and CPM transform
library(plyr)
library(dplyr)


# read in data
df_ATAC <- read.delim("/Users/wittkopp_member/Documents/ZHR_Z30_TSIM_combined_ALLclasses_ATAC_downsampled_d_merge24.txt", header = T)
df_ATAC <- df_ATAC[c((26:28),rev(2:25))]
colnames(df_ATAC) <- c("chrom", "start", "end", "ZHR_1", "ZHR_2", "ZHR_3", "Z30_1", "Z30_2", "Z30_3",
"ZHR_Z30_1_ZHR", "ZHR_Z30_2_ZHR", "ZHR_Z30_3_ZHR", "ZHR_Z30_1_Z30", "ZHR_Z30_2_Z30", "ZHR_Z30_3_Z30",
"ZHR_1_tsim", "ZHR_2_tsim", "ZHR_3_tsim", "TSIM_1", "TSIM_2", "TSIM_3",
"ZHR_TSIM_1_ZHR", "ZHR_TSIM_2_ZHR", "ZHR_TSIM_3_ZHR", "ZHR_TSIM_1_TSIM", "ZHR_TSIM_2_TSIM", "ZHR_TSIM_3_TSIM")

full_dataset <- df_ATAC

full_dataset <- full_dataset[full_dataset$ZHR_1 >= 10 & full_dataset$ZHR_2 >= 10 & full_dataset$ZHR_3 >= 10 &
full_dataset$Z30_1 >= 10 & full_dataset$Z30_2 >= 10 & full_dataset$Z30_3 >= 10 &
full_dataset$ZHR_Z30_1_ZHR >= 10 & full_dataset$ZHR_Z30_2_ZHR >= 10 & full_dataset$ZHR_Z30_3_ZHR >= 10 &
full_dataset$ZHR_Z30_1_Z30 >= 10 & full_dataset$ZHR_Z30_2_Z30 >= 10 & full_dataset$ZHR_Z30_3_Z30 >= 10 &
full_dataset$ZHR_1_tsim >= 10 & full_dataset$ZHR_2_tsim >= 10 & full_dataset$ZHR_3_tsim >= 10 &
full_dataset$TSIM_1 >= 10 & full_dataset$TSIM_2 >= 10 & full_dataset$TSIM_3 >= 10 &
full_dataset$ZHR_TSIM_1_TSIM >= 10 & full_dataset$ZHR_TSIM_2_TSIM >= 10 & full_dataset$ZHR_TSIM_3_TSIM >= 10,]

df_ATAC <- full_dataset

# CPM transformation
df_ATAC_CPM <- matrix(ncol = ncol(df_ATAC) - 3, nrow = nrow(df_ATAC)) %>% as.data.frame()

for(i in 4:ncol(df_ATAC)) {
  df_ATAC_CPM[,i] <- (df_ATAC[,i]/sum(df_ATAC[,i]))*1000000
}

df_ATAC_CPM$V1 <- df_ATAC$chrom
df_ATAC_CPM$V2 <- df_ATAC$start
df_ATAC_CPM$V3 <- df_ATAC$end
colnames(df_ATAC_CPM) <- c("chrom", "start", "end", "ZHR_1", "ZHR_2", "ZHR_3", "Z30_1", "Z30_2", "Z30_3",
"ZHR_Z30_1_ZHR", "ZHR_Z30_2_ZHR", "ZHR_Z30_3_ZHR", "ZHR_Z30_1_Z30", "ZHR_Z30_2_Z30", "ZHR_Z30_3_Z30",
"ZHR_1_tsim", "ZHR_2_tsim", "ZHR_3_tsim", "TSIM_1", "TSIM_2", "TSIM_3",
"ZHR_TSIM_1_ZHR", "ZHR_TSIM_2_ZHR", "ZHR_TSIM_3_ZHR", "ZHR_TSIM_1_TSIM", "ZHR_TSIM_2_TSIM", "ZHR_TSIM_3_TSIM")

write.table(df_ATAC_CPM, file = "/Users/wittkopp_member/Documents/ZHR_Z30_TSIM_ATAC_macs2_combined_CPM_centered_final_dm6_10min.txt", row.names = F, quote = F, sep = "\t")

# Seperate ZHR x Z30 and ZHR x TSIM to run seperately
ZHR_Z30 <- df_ATAC_CPM[c(1:15)]
ZHR_TSIM <- df_ATAC_CPM[c(1:3,16:ncol(df_ATAC_CPM))]

write.table(ZHR_Z30, file = "/Users/wittkopp_member/Documents/ZHR_Z30_sub_ATAC_macs2_combined_CPM_centered_final_dm6_20min.txt", row.names = F, quote = F, sep = "\t")
write.table(ZHR_TSIM, file = "/Users/wittkopp_member/Documents/ZHR_TSIM_sub_ATAC_macs2_combined_CPM_centered_final_dm6_20min.txt", row.names = F, quote = F, sep = "\t")
