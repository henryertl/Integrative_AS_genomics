## Upstream analyses give a matrix where:
### COL 1 = gene
### COL 2:X = ounts per constitutive exon
#### Goal: Condense exon count to a total genic count and CPM transform

library(plyr)
library(dplyr)

# read in data
df_exon_counts <- read.delim("/nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/MERGED_AND_FINAL_FILES/ATAC_RNA_integrated_files/ZHR_Z30_exon_counts_ATAC_RNA_dm6.txt", header = F)
colnames(df_exon_counts) <- c("chrom", "start", "end", "gene", "P1_1", "P1_2", "P1_3", "P1_4", "P1_5", "P1_6",
"P2_1", "P2_2", "P2_3", "P2_4", "P2_5", "P2_6",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_4_P1", "HYB_5_P1", "HYB_6_P1",
"HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "HYB_4_P2", "HYB_5_P2", "HYB_6_P2")

#re-format
df_exon_counts_gene <- df_exon_counts[,4:ncol(df_exon_counts)]

# combine exons
df_genic_counts <- ddply(df_exon_counts_gene,.(gene),summarize,P1_1=sum(P1_1),P1_2=sum(P1_2),P1_3=sum(P1_3),P1_4=sum(P1_4),P1_5=sum(P1_5),P1_6=sum(P1_6),
P2_1=sum(P2_1),P2_2=sum(P2_2),P2_3=sum(P2_3),P2_4=sum(P2_4),P2_5=sum(P2_5),P2_6=sum(P2_6),
HYB_1_P1=sum(HYB_1_P1),HYB_2_P1=sum(HYB_2_P1),HYB_3_P1=sum(HYB_3_P1),HYB_4_P1=sum(HYB_4_P1),HYB_5_P1=sum(HYB_5_P1),HYB_6_P1=sum(HYB_6_P1),
HYB_1_P2=sum(HYB_1_P2),HYB_2_P2=sum(HYB_2_P2),HYB_3_P2=sum(HYB_3_P2),HYB_4_P2=sum(HYB_4_P2),HYB_5_P2=sum(HYB_5_P2),HYB_6_P2=sum(HYB_6_P2))

full_dataset <- df_genic_counts %>% as.data.frame()

colnames(df_genic_counts) <- c("gene", "RNA_P1_1", "RNA_P1_3", "RNA_P1_5", "RNA_P2_1", "RNA_P2_3", "RNA_P2_5",
"RNA_HYB_P1_1", "RNA_HYB_P1_3", "RNA_HYB_P1_5", "RNA_HYB_P2_1", "RNA_HYB_P2_3", "RNA_HYB_P2_5",
"ATAC_P1_1", "ATAC_P1_2", "ATAC_P1_3", "ATAC_P2_1", "ATAC_P2_2", "ATAC_P2_3",
"ATAC_HYB_P1_1", "ATAC_HYB_P1_2", "ATAC_HYB_P1_3", "ATAC_HYB_P2_1", "ATAC_HYB_P2_2", "ATAC_HYB_P2_3")
write.table(df_genic_counts, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/ZHR_Z30_genic_RAWcounts_final_dm6.txt", row.names = F, quote = F, sep = "\t")

#subset >20 reads in all
full_dataset <- full_dataset[full_dataset$P1_1 >= 20 & full_dataset$P1_2 >= 20 & full_dataset$P1_3 >= 20 & full_dataset$P1_4 >= 20 & full_dataset$P1_5 >= 20 & full_dataset$P1_6 >= 20 &
full_dataset$P2_1 >= 20 & full_dataset$P2_2 >= 20 & full_dataset$P2_3 >= 20 & full_dataset$P1_4 >= 20 & full_dataset$P1_5 >= 20 & full_dataset$P1_6 >= 20 &
full_dataset$HYB_1_P1 >= 20 & full_dataset$HYB_1_P2 >= 20 & full_dataset$HYB_2_P1 >= 20 & full_dataset$HYB_2_P2 >= 20 & full_dataset$HYB_3_P1 >= 20 & full_dataset$HYB_3_P2 >=20 &
full_dataset$HYB_4_P1 >= 20 & full_dataset$HYB_4_P2 >= 20 & full_dataset$HYB_5_P1 >= 20 & full_dataset$HYB_5_P2 >= 20 & full_dataset$HYB_6_P1 >= 20 & full_dataset$HYB_6_P2 >= 20,]

#subset <1000 reads in RNA
full_dataset <- full_dataset[full_dataset$P1_1 <= 1000 & full_dataset$P1_2 <= 1000 & full_dataset$P1_3 <= 1000 & full_dataset$P1_4 <= 1000 & full_dataset$P1_5 <= 1000 & full_dataset$P1_6 <= 1000 &
full_dataset$P2_1 <= 1000 & full_dataset$P2_2 <= 1000 & full_dataset$P2_3 <= 1000 & full_dataset$P1_4 <= 1000 & full_dataset$P1_5 <= 1000 & full_dataset$P1_6 <= 1000,]

df_genic_counts <- full_dataset %>% as.data.frame()

# CPM transformation
df_genic_CPM <- matrix(ncol = ncol(df_genic_counts) - 1, nrow = nrow(df_genic_counts)) %>% as.data.frame()

for(i in 2:ncol(df_genic_counts)) {
  df_genic_CPM[,i] <- (df_genic_counts[,i]/sum(df_genic_counts[,i]))*1000000
}

df_genic_CPM$V1 <- df_genic_counts$gene
colnames(df_genic_CPM) <- c("gene", "RNA_P1_1", "RNA_P1_3", "RNA_P1_5", "RNA_P2_1", "RNA_P2_3", "RNA_P2_5",
"RNA_HYB_P1_1", "RNA_HYB_P1_3", "RNA_HYB_P1_5", "RNA_HYB_P2_1", "RNA_HYB_P2_3", "RNA_HYB_P2_5",
"ATAC_P1_1", "ATAC_P1_2", "ATAC_P1_3", "ATAC_P2_1", "ATAC_P2_2", "ATAC_P2_3",
"ATAC_HYB_P1_1", "ATAC_HYB_P1_2", "ATAC_HYB_P1_3", "ATAC_HYB_P2_1", "ATAC_HYB_P2_2", "ATAC_HYB_P2_3")
write.table(df_genic_CPM, file = "/nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/MERGED_AND_FINAL_FILES/ATAC_RNA_integrated_files/ZHR_Z30_integrate_ATACandRNA_genic_counts_CPM_final_dm6_20min.txt", row.names = F, quote = F, sep = "\t")

# split up ATAC and RNA and write respective data table files
RNA_output <- df_genic_CPM[c(1:13)]
ATAC_output <- df_genic_CPM[c(1,14:ncol(df_genic_CPM))]
ATAC_RNA_output <- join_all(ATAC_output, RNA_output, by = gene)
write.table(RNA_output, file = "/nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/MERGED_AND_FINAL_FILES/ATAC_RNA_integrated_files/ZHR_Z30_integate_RNA_genic_counts_CPM_final_dm6_20min_1000max.txt", row.names = F, quote = F, sep = "\t")
write.table(ATAC_output, file = "/nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/MERGED_AND_FINAL_FILES/ATAC_RNA_integrated_files/ZHR_Z30_integrate_ATAC_genic_counts_CPM_final_dm6_20min.txt", row.names = F, quote = F, sep = "\t")
