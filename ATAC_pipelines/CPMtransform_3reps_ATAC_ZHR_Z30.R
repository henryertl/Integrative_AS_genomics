#### Goal: Condense exon count to a total genic count and CPM transform

# read in data
df_ATAC <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Data_tables/ZHR_Z30_ATAC_merged_counts_final_macs2.bed", header = F)
colnames(df_ATAC) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3","P2_1", "P2_2", "P2_3",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")

full_dataset <- df_ATAC

full_dataset <- full_dataset[full_dataset$P1_1 >= 20 & full_dataset$P1_2 >= 20 & full_dataset$P1_3 >= 20 & full_dataset$P2_1 >= 20 & full_dataset$P2_2 >= 20 &
full_dataset$P2_3 >= 20 & full_dataset$HYB_1_P1 + full_dataset$HYB_1_P2 >= 20 & full_dataset$HYB_2_P1 + full_dataset$HYB_2_P2 >= 20 & full_dataset$HYB_3_P1 + full_dataset$HYB_3_P2 >= 20,]

df_ATAC <- full_dataset

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
write.table(df_ATAC_CPM, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Data_tables/ZHR_Z30_ATAC_macs2_CPM_final_dm6_20min.txt", row.names = F, quote = F, sep = "\t")

head(df_ATAC_CPM)

df_ATAC_CPM$diff <- df_ATAC_CPM$end - df_ATAC_CPM$start

A <- ggplot(df_ATAC_CPM, aes(x=diff)) +
geom_density()

B <- ggplot(df_ATAC_CPM, aes(x=diff)) +
geom_density()

plot_grid(A, B)
