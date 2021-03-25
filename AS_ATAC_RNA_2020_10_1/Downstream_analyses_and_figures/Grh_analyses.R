######## GOAL: join Grh motif/ChIP and variation information to Full results ###########

# load libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(rstatix)

# set wd
setwd("/Users/wittkopp_member/Code/Integrative_AS_genomics")
setwd("/Users/henryertl/Documents/Devs/Integrative_AS_genomics")

theme_main <- function() {
  theme_bw() +
  theme(
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 15),
  strip.text = element_text(size = 15),
  legend.text= element_text(size = 15),
  legend.title = element_text(size = 10),
  plot.title = element_text(size = 15, face = "bold")

)
}



# STEP 1: FIND OVERLAPPING REGIONS WITH GRH CHIP/MOTIF LOCATIONS
#### BEDTOOLS command line to prep files #####
bedtools intersect -wao -a Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.bed -b /Users/wittkopp_member/Code/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Grh_data/Grh_ChIP_motif_overlap_dm6_final.bed | uniq > Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_Grh_intersectALL.bed
bedtools intersect -wao -a Full_results_output_ZHR_TSIM_ATAC_20min_downsamp_overlap.bed -b /Users/wittkopp_member/Code/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Grh_data/Grh_ChIP_motif_overlap_dm6_final.bed > Full_results_output_ZHR_TSIM_ATAC_20min_downsamp_overlap_Grh_intersectALL.bed

## read in outputs from bedtools and clean up/add col names
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

# clean up and write table to have full results dataframe with only number of motif/ChIP occurances per site
ZHR_Z30_GRH_condensed <- ZHR_Z30_GRH[,-c(ncol(ZHR_Z30_GRH)-1, ncol(ZHR_Z30_GRH)-2, ncol(ZHR_Z30_GRH)-3, ncol(ZHR_Z30_GRH)-4, ncol(ZHR_Z30_GRH)-5)]

ZHR_Z30_GRH_condensed_count <- ddply(ZHR_Z30_GRH_condensed,.(chrom, start, end, Paste_locus, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, HYB_1_P1, HYB_2_P1, HYB_3_P1, HYB_1_P2, HYB_2_P2, HYB_3_P2,
pos_start, pos_end, P_est.mean, P_p_value,H_est.mean,H_p_value,H_P_est,H_P_p_value,P_qvalue,H_qvalue,P_H_qvalue,Direction_parent,Direction_hybrid,
P_H_ratio,Regulatory_class,trans_reg_diff,perc_cis,class,overlap_binary),nrow)

colnames(ZHR_Z30_GRH_condensed_count)[ncol(ZHR_Z30_GRH_condensed_count)] <- "count"

write.table(ZHR_Z30_GRH_condensed_count, file = "./AS_ATAC_RNA_2020_10_1/Grh_data/ZHR_Z30_GRH_condensed_count.txt", sep = "\t", quote = F, row.names = F)

# STEP 2: FIND OVERLAPPING REGIONS WITH SNP DATA IN CHIP/MOTIF LOCATIONS
#### BEDTOOLS command line to intersect with SNP data #####
bedtools intersect -wao -a Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.bed -b /Users/wittkopp_member/Code/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Grh_data/Grh_ZHR_Z30_SNPs.bed | uniq > Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_Grh_intersectALL.bed

# read in data from bedtools and add column names
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

## clean up and add column with number of SNPs in Grh ChIP/motifs per regions
ZHR_Z30_GRH_SNPs_condensed <- ZHR_Z30_GRH_SNPs[,-c(ncol(ZHR_Z30_GRH_SNPs)-1, ncol(ZHR_Z30_GRH_SNPs)-2, ncol(ZHR_Z30_GRH_SNPs)-3, ncol(ZHR_Z30_GRH_SNPs)-4, ncol(ZHR_Z30_GRH_SNPs)-5)]

ZHR_Z30_GRH_SNPs_condensed_count <- ddply(ZHR_Z30_GRH_SNPs_condensed,.(chrom, start, end, Paste_locus, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, HYB_1_P1, HYB_2_P1, HYB_3_P1, HYB_1_P2, HYB_2_P2, HYB_3_P2,
pos_start, pos_end, P_est.mean, P_p_value,H_est.mean,H_p_value,H_P_est,H_P_p_value,P_qvalue,H_qvalue,P_H_qvalue,Direction_parent,Direction_hybrid,
P_H_ratio,Regulatory_class,trans_reg_diff,perc_cis,class,overlap_SNP_binary),nrow)

colnames(ZHR_Z30_GRH_SNPs_condensed_count)[ncol(ZHR_Z30_GRH_SNPs_condensed_count)] <- "count_SNP"

## clean up and join grh motif (step 1) and SNP dfs (step 2)
ZHR_Z30_GRH_SNPs_condensed_count_mintojoin <- ZHR_Z30_GRH_SNPs_condensed_count[,c(4,35,36)]

# read in output for STEP 1
ZHR_Z30_GRH_condensed_count <- read.delim("./AS_ATAC_RNA_2020_10_1/Grh_data/ZHR_Z30_GRH_condensed_count.txt", header = T)
ZHR_Z30_GRH_SNPs_final <- join_all(list(ZHR_Z30_GRH_SNPs_condensed_count_mintojoin, ZHR_Z30_GRH_condensed_count), type = "full", by = "Paste_locus")

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

write.table(ZHR_Z30_GRH_SNPs_final, file = "./AS_ATAC_RNA_2020_10_1/Grh_data/ZHR_Z30_GRH_SNPs_count.txt", sep = "\t", quote = F, row.names = F)

# STEP 3: add total SNP info to df
##### BEDtools to get total SNPs ######
bedtools intersect -wao -a Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.bed -b /Users/wittkopp_member/Code/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/zhr_z30_SNPs_dm6_genome.bed | uniq > Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_SNPs_intersectALL.bed

# read in data from bedtools and add column names/clean up
ZHR_Z30_regions_ALL_SNPs <- read.delim("./AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_SNPs_intersectALL.bed", header = F)
colnames(ZHR_Z30_regions_ALL_SNPs) <- c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2",
"pos_start", "pos_end", "P_est.mean", "P_p_value","H_est.mean","H_p_value","H_P_est","H_P_p_value","P_qvalue","H_qvalue","P_H_qvalue","Direction_parent","Direction_hybrid",
"P_H_ratio","Regulatory_class","trans_reg_diff","perc_cis","class", "chrom_snp", "start_snp","end_snp","overlap_snp")

## clean up and add column with number of total SNPs in region
ZHR_Z30_regions_ALL_SNPs <- ZHR_Z30_regions_ALL_SNPs[,c(1:(ncol(ZHR_Z30_regions_ALL_SNPs) - 4))]

ZHR_Z30_regions_ALL_SNPs <- ddply(ZHR_Z30_regions_ALL_SNPs,.(chrom, start, end, Paste_locus, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, HYB_1_P1, HYB_2_P1, HYB_3_P1, HYB_1_P2, HYB_2_P2, HYB_3_P2,
pos_start, pos_end, P_est.mean, P_p_value,H_est.mean,H_p_value,H_P_est,H_P_p_value,P_qvalue,H_qvalue,P_H_qvalue,Direction_parent,Direction_hybrid,
P_H_ratio,Regulatory_class,trans_reg_diff,perc_cis,class),nrow)

ZHR_Z30_regions_ALL_SNPs_final <- ZHR_Z30_regions_ALL_SNPs[,c(4,ncol(ZHR_Z30_regions_ALL_SNPs))]
colnames(ZHR_Z30_regions_ALL_SNPs_final) <- c("Paste_locus", "num_SNPs_total")

# join with Grh_snp file
ZHR_Z30_regions_ALL_SNPs_GRh_final <- join_all(list(ZHR_Z30_regions_ALL_SNPs_final, ZHR_Z30_GRH_SNPs_final), type = 'full', by = 'Paste_locus') %>% unique() %>% na.omit() %>% as.data.frame()
ZHR_Z30_regions_ALL_SNPs_GRh_final$SNPs_motif_norm <- NULL

colnames(ZHR_Z30_regions_ALL_SNPs_GRh_final) <- c("Paste_locus", "total_region_snps", "overlap_grh_snp_binary", "grh_snp_count", "chrom", "start", "end", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2",
"pos_start", "pos_end", "P_est.mean", "P_p_value","H_est.mean","H_p_value","H_P_est","H_P_p_value","P_qvalue","H_qvalue","P_H_qvalue","Direction_parent","Direction_hybrid",
"P_H_ratio","Regulatory_class","trans_reg_diff","perc_cis","class", "overlap_grh_motif_binary", "grh_motif_counts")

ZHR_Z30_regions_ALL_SNPs_GRh_final$P_est.mean_abs <- abs(ZHR_Z30_regions_ALL_SNPs_GRh_final$P_est.mean)

write.table(ZHR_Z30_regions_ALL_SNPs_GRh_final, file = "./AS_ATAC_RNA_2020_10_1/Grh_data/ZHR_Z30_GRH_SNPs_count_FINALALL.txt", sep = "\t", quote = F, row.names = F)

ZHR_Z30_regions_ALL_SNPs_GRh_final <- read.delim("./AS_ATAC_RNA_2020_10_1/Grh_data/ZHR_Z30_GRH_SNPs_count_FINALALL.txt", header = T)

# STEP 4: ANALYSES
ZHR_Z30_regions_ALL_SNPs_GRh_final$class <- factor(ZHR_Z30_regions_ALL_SNPs_GRh_final$class,levels = c("start", "end", "inter", "intra"))

ZHR_Z30_regions_ALL_SNPs_GRh_final_grh_overlap <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP",]
ZHR_Z30_regions_ALL_SNPs_GRh_final_sig_sub <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$P_qvalue < 0.05 | ZHR_Z30_regions_ALL_SNPs_GRh_final$H_qvalue < 0.05,]
ZHR_Z30_regions_ALL_SNPs_GRh_final_grh_overlap_7_9SNPs <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP" & ZHR_Z30_regions_ALL_SNPs_GRh_final$total_region_snps > 4 & ZHR_Z30_regions_ALL_SNPs_GRh_final$total_region_snps < 12,]

# are regions with grh motifs more or less divergent?
G <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$P_qvalue < 0.05 | ZHR_Z30_regions_ALL_SNPs_GRh_final$H_qvalue < 0.05,] %>%
ggplot(aes(x=overlap_grh_motif_binary, y=abs(P_est.mean), fill=overlap_grh_motif_binary))+
geom_boxplot(notch=T) +
ylim(0,1)

pairwise.wilcox.test(ZHR_Z30_regions_ALL_SNPs_GRh_final$P_est.mean_abs,ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary)
## regions are not more or less divergent  - q: 0.44

# are regions with ChIP grh binding and motifs different from those without?
H <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$P_qvalue < 0.05 | ZHR_Z30_regions_ALL_SNPs_GRh_final$H_qvalue < 0.05,] %>%
ggplot(aes(x=overlap_grh_motif_binary, y=perc_cis, fill=overlap_grh_motif_binary))+
geom_boxplot(notch = T)

pairwise.wilcox.test(ZHR_Z30_regions_ALL_SNPs_GRh_final_sig_sub$perc_cis,ZHR_Z30_regions_ALL_SNPs_GRh_final_sig_sub$overlap_grh_motif_binary)
## percent cis is slightly higher for regions with grh binding - q: 0.012

# is there a relationship between number of grh motifs and percent cis?
ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_counts > 0,] %>%
ggplot() +
geom_density(aes(x=grh_motif_counts))

ZHR_Z30_regions_ALL_SNPs_GRh_final <- as.data.frame(ZHR_Z30_regions_ALL_SNPs_GRh_final)

ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_bin <- "NA"

for (i in 1:nrow(ZHR_Z30_regions_ALL_SNPs_GRh_final)) {

if (ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_counts[i] < 1){

	ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_bin[i] <- "0"

} else if (ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_counts[i] > 0 & ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_counts[i] < 3){

	ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_bin[i] <- "1-2"

} else if (ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_counts[i] > 2 & ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_counts[i] < 8){

	ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_bin[i] <- "3-7"

} else if (ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_counts[i] > 7){

	ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_bin[i] <- ">7"

}
}

nrow(ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_bin == ">7",])

ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_bin <- factor(ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_bin,levels = c("0", "1-2", "3-7", ">7"))

G <- ZHR_Z30_regions_ALL_SNPs_GRh_final %>%
ggplot(aes(x=grh_motif_bin, y=P_est.mean_abs, fill=grh_motif_bin))+
geom_boxplot(notch=T) +
ylim(0,2)

H <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$P_qvalue < 0.05 | ZHR_Z30_regions_ALL_SNPs_GRh_final$H_qvalue < 0.05,] %>%
ggplot(aes(x=grh_motif_bin, y=perc_cis, fill=grh_motif_bin))+
geom_boxplot(notch=T)

I <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$P_qvalue < 0.05 | ZHR_Z30_regions_ALL_SNPs_GRh_final$H_qvalue < 0.05,] %>%
ggplot(aes(x=grh_motif_bin, y=total_region_snps, fill=grh_motif_bin))+
geom_boxplot(notch=T)

I <- ZHR_Z30_regions_ALL_SNPs_GRh_final %>%
ggplot(aes(x=class, y=total_region_snps, fill=class))+
geom_boxplot(notch=TRUE) +
theme_main() +
ylab("Number of SNPs in region") +
xlab("") +
scale_fill_discrete(guide=FALSE) +
scale_x_discrete(labels=c("txStart","txEnd", "intergenic", "intragenic")) +
ylim(0,30)
ggsave(I, file = "./AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/ZHR_Z30_Full_results_output_ALL_classes_total_SNPs.pdf")


ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_snps_per_grh_motif_norm_by_total_snps <- (ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_snp_count / ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_counts)/ZHR_Z30_regions_ALL_SNPs_GRh_final$total_region_snps

J <- ZHR_Z30_regions_ALL_SNPs_GRh_final %>%
ggplot(aes(x=grh_motif_bin, y=log(grh_snps_per_grh_motif_norm_by_total_snps), fill=grh_motif_bin))+
geom_boxplot(notch = T)

####
ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$P_est.mean_abs > 0.5,] %>%
ggplot() +
geom_point(aes(x=perc_cis, y=abs(P_est.mean)))

G <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$P_qvalue < 0.05 | ZHR_Z30_regions_ALL_SNPs_GRh_final$H_qvalue < 0.05 & ZHR_Z30_regions_ALL_SNPs_GRh_final$class == "inter",] %>%
ggplot(aes(x=grh_motif_bin, y=perc_cis, fill=grh_motif_bin))+
geom_boxplot(notch=T) +
ylim(0,1)

pairwise.wilcox.test(ZHR_Z30_regions_ALL_SNPs_GRh_final$perc_cis,ZHR_Z30_regions_ALL_SNPs_GRh_final$grh_motif_bin)


# do variable grh motifs have more divergence than nonvariable?
B <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP",] %>%
ggplot(aes(x=overlap_grh_snp_binary, y=abs(P_est.mean), fill=overlap_grh_snp_binary))+
geom_boxplot(notch=T) +
ylim(0,1)

pairwise.wilcox.test(ZHR_Z30_regions_ALL_SNPs_GRh_final_grh_overlap$P_est.mean_abs,ZHR_Z30_regions_ALL_SNPs_GRh_final_grh_overlap$overlap_grh_snp_binary)
## yes however this could be due to variable grh motifs being in overall more variable regions - q: 0.018

# are variable grh motifs in overall more variable regions?
C <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP",] %>%
ggplot(aes(x=overlap_grh_snp_binary, y=total_region_snps, fill=overlap_grh_snp_binary))+
geom_boxplot(notch=F)

pairwise.wilcox.test(ZHR_Z30_regions_ALL_SNPs_GRh_final_grh_overlap$total_region_snps,ZHR_Z30_regions_ALL_SNPs_GRh_final_grh_overlap$overlap_grh_snp_binary)
## yes these regions are just more variable... need to isolate the input from grh variation - q: <2e-16


# compare divergence between grh variation/no variation for regions with equalized total varation
## first look at distr of total region snps
D <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP",] %>%
ggplot(aes(x=total_region_snps))+
geom_density()
## mean is ~8

E <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP",] %>%
ggplot(aes(x=overlap_grh_snp_binary, y=abs(P_est.mean), fill=overlap_grh_snp_binary))+
geom_boxplot(notch=T)
pairwise.wilcox.test(ZHR_Z30_regions_ALL_SNPs_GRh_final_grh_overlap$P_est.mean_abs,ZHR_Z30_regions_ALL_SNPs_GRh_final_grh_overlap$overlap_grh_snp_binary)

## if only subset regions w/ 8 total snps, then there is a slight difference for regions with variable grh motifs

G <- ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$overlap_grh_motif_binary == "OVERLAP" & ZHR_Z30_regions_ALL_SNPs_GRh_final$total_region_snps == 8,] %>%
ggplot(aes(x=overlap_grh_snp_binary, y=perc_cis, fill=overlap_grh_snp_binary))+
geom_boxplot(notch=T)
## pattern above is exxgerated if extend range of permitted total snp counts
