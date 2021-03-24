
library(plyr)
library(ggplot2)

# set wd
setwd("/Users/wittkopp_member/Code/Integrative_AS_genomics")

#### BEDTOOLS command line to prep files #####
bedtools intersect -wao -a Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.bed -b /Users/wittkopp_member/Code/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Grh_data/Grh_ChIP_motif_overlap_dm6_final.bed | uniq > Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_Grh_intersectALL.bed
bedtools intersect -wao -a Full_results_output_ZHR_TSIM_ATAC_20min_downsamp_overlap.bed -b /Users/wittkopp_member/Code/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Grh_data/Grh_ChIP_motif_overlap_dm6_final.bed > Full_results_output_ZHR_TSIM_ATAC_20min_downsamp_overlap_Grh_intersectALL.bed

# read in outputs and clean up/add col names
ZHR_Z30_GRH <- read.delim("./AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_Grh_intersectALL.bed", header = F) %>% as.data.frame()
colnames(ZHR_Z30_GRH) <- c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2",
"pos_start", "pos_end", "P_est.mean", "P_p_value","H_est.mean","H_p_value","H_P_est","H_P_p_value","P_qvalue","H_qvalue","P_H_qvalue","Direction_parent","Direction_hybrid",
"P_H_ratio","Regulatory_class","trans_reg_diff","perc_cis","class", "chrom_grh", "start_grh","end_grh","score_grh","overlap")

# loop to reassign overlap as YES or NO
ZHR_Z30_GRH$overlap_binary <- "NA"

for (i in 1:nrow(ZHR_Z30_GRH)) {

if (ZHR_Z30_GRH$overlap[i] != 0){

	ZHR_Z30_GRH$overlap_binary[i] <- "OVERLAP"

} else if (ZHR_Z30_GRH$overlap[i] == 0){

	ZHR_Z30_GRH$overlap_binary[i] <- "NO_OVERLAP"

}
}

ZHR_Z30_GRH_condensed <- ZHR_Z30_GRH[,-c(ncol(ZHR_Z30_GRH)-1, ncol(ZHR_Z30_GRH)-2, ncol(ZHR_Z30_GRH)-3, ncol(ZHR_Z30_GRH)-4, ncol(ZHR_Z30_GRH)-5)]

ZHR_Z30_GRH_condensed_count <- ddply(ZHR_Z30_GRH_condensed,.(chrom, start, end, Paste_locus, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, HYB_1_P1, HYB_2_P1, HYB_3_P1, HYB_1_P2, HYB_2_P2, HYB_3_P2,
pos_start, pos_end, P_est.mean, P_p_value,H_est.mean,H_p_value,H_P_est,H_P_p_value,P_qvalue,H_qvalue,P_H_qvalue,Direction_parent,Direction_hybrid,
P_H_ratio,Regulatory_class,trans_reg_diff,perc_cis,class,overlap_binary),nrow)

colnames(ZHR_Z30_GRH_condensed_count)[ncol(ZHR_Z30_GRH_condensed_count)] <- "count"

write.table(ZHR_Z30_GRH_condensed_count, file = "./AS_ATAC_RNA_2020_10_1/Grh_data/ZHR_Z30_GRH_condensed_count.txt", sep = "\t", quote = F, row.names = F)

#### BEDTOOLS command line to intersect with SNP data #####
bedtools intersect -wao -a Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.bed -b /Users/wittkopp_member/Code/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Grh_data/Grh_ZHR_Z30_SNPs.bed | uniq > Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_Grh_intersectALL.bed
####

ZHR_Z30_GRH_SNPs <- read.delim("./AS_ATAC_RNA_2020_10_1/Grh_data/ZHR_Z30_GRH_condensed_count_SNPintersect.bed", header = F)

colnames(ZHR_Z30_GRH_SNPs) <- c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2",
"pos_start", "pos_end", "P_est.mean", "P_p_value","H_est.mean","H_p_value","H_P_est","H_P_p_value","P_qvalue","H_qvalue","P_H_qvalue","Direction_parent","Direction_hybrid",
"P_H_ratio","Regulatory_class","trans_reg_diff","perc_cis","class", "chrom_grh_snp", "start_grh_snp","end_grh_snp","score_grh_snp","overlap_SNP")

# loop to reassign overlap as YES or NO
ZHR_Z30_GRH_SNPs$overlap_SNP_binary <- "NA"

for (i in 1:nrow(ZHR_Z30_GRH_SNPs)) {

if (ZHR_Z30_GRH_SNPs$overlap_SNP[i] != 0){

	ZHR_Z30_GRH_SNPs$overlap_SNP_binary[i] <- "OVERLAP"

} else if (ZHR_Z30_GRH_SNPs$overlap_SNP[i] == 0){

	ZHR_Z30_GRH_SNPs$overlap_SNP_binary[i] <- "NO_OVERLAP"

}
}

ZHR_Z30_GRH_SNPs_condensed <- ZHR_Z30_GRH_SNPs[,-c(ncol(ZHR_Z30_GRH_SNPs)-1, ncol(ZHR_Z30_GRH_SNPs)-2, ncol(ZHR_Z30_GRH_SNPs)-3, ncol(ZHR_Z30_GRH_SNPs)-4, ncol(ZHR_Z30_GRH_SNPs)-5)]

ZHR_Z30_GRH_SNPs_condensed_count <- ddply(ZHR_Z30_GRH_SNPs_condensed,.(chrom, start, end, Paste_locus, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, HYB_1_P1, HYB_2_P1, HYB_3_P1, HYB_1_P2, HYB_2_P2, HYB_3_P2,
pos_start, pos_end, P_est.mean, P_p_value,H_est.mean,H_p_value,H_P_est,H_P_p_value,P_qvalue,H_qvalue,P_H_qvalue,Direction_parent,Direction_hybrid,
P_H_ratio,Regulatory_class,trans_reg_diff,perc_cis,class,overlap_SNP_binary),nrow)

colnames(ZHR_Z30_GRH_SNPs_condensed_count)[ncol(ZHR_Z30_GRH_SNPs_condensed_count)] <- "count_SNP"

# clean up and join grh motif and SNP dfs

ZHR_Z30_GRH_SNPs_condensed_count_mintojoin <- ZHR_Z30_GRH_SNPs_condensed_count[,c(4,35,36)]

ZHR_Z30_GRH_SNPs_final <- join_all(list(ZHR_Z30_GRH_SNPs_condensed_count_mintojoin, ZHR_Z30_GRH_condensed_count), type = "full", by = "Paste_locus")

ZHR_Z30_GRH_SNPs_final$SNPs_motif_norm <- ZHR_Z30_GRH_SNPs_final$count_SNP / ZHR_Z30_GRH_SNPs_final$count

# loops to clean up and add 0s where needed

for (i in 1:nrow(ZHR_Z30_GRH_SNPs_final)) {

if (ZHR_Z30_GRH_SNPs_final$overlap_binary[i] == "NO_OVERLAP"){

	ZHR_Z30_GRH_SNPs_final$count[i] <- 0

}
}

for (i in 1:nrow(ZHR_Z30_GRH_SNPs_final)) {

if (ZHR_Z30_GRH_SNPs_final$overlap_SNP_binary[i] == "NO_OVERLAP"){

	ZHR_Z30_GRH_SNPs_final$count_SNP[i] <- 0

}
}

for (i in 1:nrow(ZHR_Z30_GRH_SNPs_final)) {

if (ZHR_Z30_GRH_SNPs_final$overlap_binary[i] == "NO_OVERLAP"){

	ZHR_Z30_GRH_SNPs_final$SNPs_motif_norm[i] <- 0

}
}

for (i in 1:nrow(ZHR_Z30_GRH_SNPs_final)) {

if (ZHR_Z30_GRH_SNPs_final$overlap_SNP_binary[i] == "NO_OVERLAP"){

	ZHR_Z30_GRH_SNPs_final$SNPs_motif_norm[i] <- 0

}
}



A <- ZHR_Z30_GRH_SNPs_final[ZHR_Z30_GRH_SNPs_final$overlap_binary == "OVERLAP",] %>%
ggplot(aes(x=overlap_SNP_binary, y=perc_cis, fill=overlap_SNP_binary)) +
geom_boxplot(notch=T)



##### get total SNPs in regions to ensure no difference in SNPs ######
bedtools intersect -wao -a Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.bed -b /Users/wittkopp_member/Code/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/zhr_z30_SNPs_dm6_genome.bed | uniq > Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_SNPs_intersectALL.bed


ZHR_Z30_regions_ALL_SNPs <- read.delim("./AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_SNPs_intersectALL.bed", header = F)
colnames(ZHR_Z30_regions_ALL_SNPs) <- c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2",
"pos_start", "pos_end", "P_est.mean", "P_p_value","H_est.mean","H_p_value","H_P_est","H_P_p_value","P_qvalue","H_qvalue","P_H_qvalue","Direction_parent","Direction_hybrid",
"P_H_ratio","Regulatory_class","trans_reg_diff","perc_cis","class", "chrom_snp", "start_snp","end_snp","overlap_snp")


ZHR_Z30_regions_ALL_SNPs <- ZHR_Z30_regions_ALL_SNPs[,c(1:(ncol(ZHR_Z30_regions_ALL_SNPs) - 4))]

ZHR_Z30_regions_ALL_SNPs <- ddply(ZHR_Z30_regions_ALL_SNPs,.(chrom, start, end, Paste_locus, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, HYB_1_P1, HYB_2_P1, HYB_3_P1, HYB_1_P2, HYB_2_P2, HYB_3_P2,
pos_start, pos_end, P_est.mean, P_p_value,H_est.mean,H_p_value,H_P_est,H_P_p_value,P_qvalue,H_qvalue,P_H_qvalue,Direction_parent,Direction_hybrid,
P_H_ratio,Regulatory_class,trans_reg_diff,perc_cis,class),nrow)

ZHR_Z30_regions_ALL_SNPs_final <- ZHR_Z30_regions_ALL_SNPs[,c(4,ncol(ZHR_Z30_regions_ALL_SNPs))]
colnames(ZHR_Z30_regions_ALL_SNPs_final) <- c("Paste_locus", "num_SNPs_total")

# join with Grh_snp file
ZHR_Z30_regions_ALL_SNPs_GRh_final <- join_all(list(ZHR_Z30_regions_ALL_SNPs_final, ZHR_Z30_GRH_SNPs_final), type = 'full', by = 'Paste_locus') %>% unique() %>% na.omit() %>% as.data.frame()

ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_snps_motif_totalsnp_norm <- ZHR_Z30_regions_ALL_SNPs_GRh_final$SNPs_motif_norm / ZHR_Z30_regions_ALL_SNPs_GRh_final$num_SNPs_total
colnames(ZHR_Z30_regions_ALL_SNPs_GRh_final) <- c("Paste_locus", "total_region_snps", "overlap_grh_snp_binary", "grh_snp_count", "chrom", "start", "end", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2",
"pos_start", "pos_end", "P_est.mean", "P_p_value","H_est.mean","H_p_value","H_P_est","H_P_p_value","P_qvalue","H_qvalue","P_H_qvalue","Direction_parent","Direction_hybrid",
"P_H_ratio","Regulatory_class","trans_reg_diff","perc_cis","class", "overlap_grh_motif_binary", "grh_motif_counts", "grh_snps_norm_by_motif_count", "grh_snps_norm_by_motif_count_total_SNPs")

# is there a relationship between grh motif presence and divergence?

A <- ZHR_Z30_regions_ALL_SNPs_GRh_final %>%
ggplot(aes(x=overlap_grh_motif_binary, y=abs(P_est.mean), fill=overlap_grh_motif_binary))+
geom_boxplot(notch=F) +
ylim(0,1)
## not really

# do variable grh motifs have more divergence than nonvariable?
B <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP",] %>%
ggplot(aes(x=overlap_grh_snp_binary, y=abs(P_est.mean), fill=overlap_grh_snp_binary))+
geom_boxplot(notch=F) +
ylim(0,1)
## yes however this could be due to variable grh motifs being in overall more variable regions

# are variable grh motifs in overall more variable regions?
C <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP",] %>%
ggplot(aes(x=overlap_grh_snp_binary, y=total_region_snps, fill=overlap_grh_snp_binary))+
geom_boxplot(notch=F)
## yes these regions are just more variable... need to isolate the input from grh variation

# compare divergence between grh variation/no variation for regions with equalized total varation
## first look at distr of total region snps
D <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP",] %>%
ggplot(aes(x=total_region_snps))+
geom_density()
## mean is ~8

E <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP" & ZHR_Z30_regions_ALL_SNPs_GRh_final$total_region_snps == 8,] %>%
ggplot(aes(x=overlap_grh_snp_binary, y=abs(P_est.mean), fill=overlap_grh_snp_binary))+
geom_boxplot(notch=F)
## if only subset regions w/ 8 total snps, then there is a slight difference for regions with variable grh motifs

F <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP" & ZHR_Z30_regions_ALL_SNPs_GRh_final$total_region_snps > 6 & ZHR_Z30_regions_ALL_SNPs_GRh_final$total_region_snps < 12,] %>%
ggplot(aes(x=overlap_grh_snp_binary, y=perc_cis, fill=overlap_grh_snp_binary))+
geom_boxplot()
## pattern above is exxgerated if extend range of permitted total snp counts

# are regions with ChIP grh binding and motifs different from those without?
G <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$P_qvalue < 0.05,] %>%
ggplot(aes(x=overlap_grh_motif_binary, y=perc_cis, fill=overlap_grh_motif_binary))+
geom_boxplot()
## percent cis is slightly higher for regions with grh binding

H <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$P_qvalue < 0.05,] %>%
ggplot(aes(x=overlap_grh_motif_binary, y=abs(P_est.mean), fill=overlap_grh_motif_binary))+
geom_boxplot()
## regions are slightly less divergent with grh binding
