## Upstream analyses give a matrix where:
### COL 1 = gene
### COL 2:X = gene counts per constitutive exon
#### Goal: Condense exon count to a total genic count and CPM transform

library(dplyr)

# read in data
df_exon_counts <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/ZHR_Z30_genic_counts_dm6.txt", header = F)
colnames(df_exon_counts) <- c("chrom", "start", "end", "gene", "P1_1", "P1_2", "P1_3", "P1_4", "P1_5", "P1_6",
"P2_1", "P2_2", "P2_3", "P2_4", "P2_5", "P2_6",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_4_P1", "HYB_5_P1", "HYB_6_P1",
"HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "HYB_4_P2", "HYB_5_P2", "HYB_6_P2")

#re-format
df_exon_counts_gene <- df_exon_counts[,4:ncol(df_exon_counts)]

# combine exons
df_RNA_data_genic_counts <- ddply(df_exon_counts_gene,.(gene),summarize,P1_1=sum(P1_1),P1_2=sum(P1_2),P1_3=sum(P1_3),P1_4=sum(P1_4),P1_5=sum(P1_5),P1_6=sum(P1_6),
P2_1=sum(P2_1),P2_2=sum(P2_2),P2_3=sum(P2_3),P2_4=sum(P2_4),P2_5=sum(P2_5),P2_6=sum(P2_6),
HYB_1_P1=sum(HYB_1_P1),HYB_2_P1=sum(HYB_2_P1),HYB_3_P1=sum(HYB_3_P1),HYB_4_P1=sum(HYB_4_P1),HYB_5_P1=sum(HYB_5_P1),HYB_6_P1=sum(HYB_6_P1),
HYB_1_P2=sum(HYB_1_P2),HYB_2_P2=sum(HYB_2_P2),HYB_3_P2=sum(HYB_3_P2),HYB_4_P2=sum(HYB_4_P2),HYB_5_P2=sum(HYB_5_P2),HYB_6_P2=sum(HYB_6_P2))
write.table(df_RNA_data_genic_counts, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/ZHR_Z30_genic_counts_final_dm6.txt", row.names = F, quote = F, sep = "\t")

#subset >20 reads and <1000 in all
full_dataset <- df_RNA_data_genic_counts %>% as.data.frame()

full_dataset <- full_dataset[full_dataset$P1_1 >= 20 & full_dataset$P1_2 >= 20 & full_dataset$P1_3 >= 20 & full_dataset$P1_4 >= 20 & full_dataset$P1_5 >= 20 & full_dataset$P1_6 >= 20 &
full_dataset$P2_1 >= 20 & full_dataset$P2_2 >= 20 & full_dataset$P2_3 >= 20 & full_dataset$P1_4 >= 20 & full_dataset$P1_5 >= 20 & full_dataset$P1_6 >= 20 &
full_dataset$HYB_1_P1 >= 20 & full_dataset$HYB_1_P2 >= 20 & full_dataset$HYB_2_P1 >= 20 & full_dataset$HYB_2_P2 >= 20 & full_dataset$HYB_3_P1 >= 20 & full_dataset$HYB_3_P2 >=20 &
full_dataset$HYB_4_P1 >= 20 & full_dataset$HYB_4_P2 >= 20 & full_dataset$HYB_5_P1 >= 20 & full_dataset$HYB_5_P2 >= 20 & full_dataset$HYB_6_P1 >= 20 & full_dataset$HYB_6_P2 >= 20,]

full_dataset <- full_dataset[full_dataset$P1_1 <= 1000 & full_dataset$P1_2 <= 1000 & full_dataset$P1_3 <= 1000 & full_dataset$P1_4 <= 1000 & full_dataset$P1_5 <= 1000 & full_dataset$P1_6 <= 1000 &
full_dataset$P2_1 <= 1000 & full_dataset$P2_2 <= 1000 & full_dataset$P2_3 <= 1000 & full_dataset$P1_4 <= 1000 & full_dataset$P1_5 <= 1000 & full_dataset$P1_6 <= 1000 &
full_dataset$HYB_1_P1 <= 1000 & full_dataset$HYB_1_P2 <= 1000 & full_dataset$HYB_2_P1 <= 1000 & full_dataset$HYB_2_P2 <= 1000 & full_dataset$HYB_3_P1 <= 1000 & full_dataset$HYB_3_P2 <=1000 &
full_dataset$HYB_4_P1 <= 1000 & full_dataset$HYB_4_P2 <= 1000 & full_dataset$HYB_5_P1 <= 1000 & full_dataset$HYB_5_P2 <= 1000 & full_dataset$HYB_6_P1 <= 1000 & full_dataset$HYB_6_P2 <= 1000,]

df_RNA_data_genic_counts <- full_dataset %>% as.data.frame()

# CPM transformation
df_RNA_data_genic_CPM <- matrix(ncol = ncol(df_RNA_data_genic_counts) - 1, nrow = nrow(df_RNA_data_genic_counts)) %>% as.data.frame()

for(i in 2:ncol(df_RNA_data_genic_counts)) {
  df_RNA_data_genic_CPM[,i] <- (df_RNA_data_genic_counts[,i]/sum(df_RNA_data_genic_counts[,i]))*1000000
}

df_RNA_data_genic_CPM$V1 <- df_RNA_data_genic_counts$gene
colnames(df_RNA_data_genic_CPM) <- c("gene", "P1_1", "P1_2", "P1_3", "P1_4", "P1_5", "P1_6",
"P2_1", "P2_2", "P2_3", "P2_4", "P2_5", "P2_6",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_4_P1", "HYB_5_P1", "HYB_6_P1",
"HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "HYB_4_P2", "HYB_5_P2", "HYB_6_P2")
write.table(df_RNA_data_genic_CPM, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/ZHR_Z30_genic_counts_CPM_final_dm6_20min_1000max.txt", row.names = F, quote = F, sep = "\t")


df_RNA_data_genic_CPM$test <- log2(df_RNA_data_genic_CPM$P1_6/df_RNA_data_genic_CPM$P2_6)

ggplot(full_dataset, aes(x=P1_1, y=P2_1)) +
geom_point()
