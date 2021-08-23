# load libraries and functions
library(stringr)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(rstatix)

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

rm(list = ls())

setwd("/Users/henryertl/Documents/Devs")

# Read in Bayes output file
ALL <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.txt", header = T)

#######################################################################
################### Add on grainyhead annotations #####################
#######################################################################
## get diff accessible regions based on grh activity

### run diffbind
library(DiffBind)
setwd("/Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles")
df <- read.csv("./Grh_KO_cntrl_diffbind.csv", header = TRUE)
df_dbind <- dba(sampleSheet = df)
df_dbind <- dba.count(df_dbind)
df_dbind <- dba.analyze(df_dbind)
df_dbind.DB <- dba.report(df_dbind)
write.table(df_dbind.DB, file = "./grh_KO_cntrl_diffbind.txt", quote = F, row.names = F, sep = "\t")

### intersect with my ATAC regions
bedtools intersect -wao -a /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.bed -b /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/grh_KO_cntrl_diffbind.bed | uniq > /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.bed

### read in intersected output from above
ZHR_Z30_GRH_diffbind <- read.delim("/Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Grh_data/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_Grh_intersectALL_diffbind.bed", header = F) %>% as.data.frame()
colnames(ZHR_Z30_GRH_diffbind) <- c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2",
"pos_start", "pos_end", "P_est.mean", "P_p_value","H_est.mean","H_p_value","H_P_est","H_P_p_value","P_qvalue","H_qvalue","P_H_qvalue","Direction_parent","Direction_hybrid",
"P_H_ratio","Regulatory_class","trans_reg_diff","perc_cis","class", "chrom_grh", "start_grh","end_grh","overlap")
# loop to reassign overlap as YES or NO
ZHR_Z30_GRH_diffbind$overlap_binary <- "NA"
for (i in 1:nrow(ZHR_Z30_GRH_diffbind)) {
if (ZHR_Z30_GRH_diffbind$overlap[i] != 0){
	ZHR_Z30_GRH_diffbind$overlap_binary[i] <- "OVERLAP"
} else if (ZHR_Z30_GRH_diffbind$overlap[i] == 0){
	ZHR_Z30_GRH_diffbind$overlap_binary[i] <- "NO_OVERLAP"
}
}
# cleanup
ZHR_Z30_GRH_diffbind_tojoin <- ZHR_Z30_GRH_diffbind[,c(4,ncol(ZHR_Z30_GRH_diffbind))]

# join to main df
ALL <- left_join(ALL, ZHR_Z30_GRH_diffbind_tojoin, by = "Paste_locus") %>% unique() %>% na.omit()

#######################################################################
################### Add on promoter annotations #####################
#######################################################################

# intersect my ATAC regions with downloaded annotations
bedtools intersect -wao -a /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.txt -b /Users/henryertl/Downloads/Promoters_characterization.bed > Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_promoters.bed

# clean up and join to full df
df <- read.delim("/Users/henryertl/Downloads/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_promoters.bed", header = F)
df <- df[df$V50 != 0,]
df <- df[,c(4,40)] %>% unique()
colnames(df) <- c("Paste_locus", "Promoter_type")
ALL <- join_all(list(ALL, df), by = "Paste_locus") %>% unique()
ALL$Promoter_type[is.na(ALL$Promoter_type)] <- "U"

## factor sort to have classes in same order in every plot
ALL$class <- factor(ALL$class,levels = c("start", "end", "inter", "intra"))

write.table(ALL, file = "/Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/ZHR_Z30_ATAC_Full_results_all_annotations.txt", row.names=F, quote= F, sep="\t")

# Plot % cis across species
### perc cis ####
R <- ALL[ALL$P_qvalue < 0.05 | ALL$H_qvalue < 0.05,] %>%
ggplot(aes(x=class, y=perc_cis, fill=class)) +
geom_boxplot(notch=TRUE) +
theme_main() +
ylab("Percent cis") +
xlab("") +
scale_fill_discrete(guide=FALSE) +
scale_x_discrete(labels=c("txStart","txEnd", "intergenic", "intragenic"))
ggsave(R, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/ZHR_Z30_Full_results_output_ALL_classes_perc_cis.pdf")

pairwise.wilcox.test(ALL$perc_cis,ALL$class)

### magnitude of divergence for differentially expressed ####
S <- ALL[ALL$P_qvalue < 0.05 | ALL$H_qvalue < 0.05,] %>%
ggplot(aes(x=class, y=abs(P_est.mean), fill=class)) +
geom_boxplot(notch=TRUE) +
theme_main() +
ylim(0,1) +
ylab("Estimated accessibility divergence") +
xlab("") +
scale_fill_discrete(guide=FALSE) +
scale_x_discrete(labels=c("txStart","txEnd", "intergenic", "intragenic"))
ggsave(S, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/ZHR_Z30_Full_results_output_ALL_classes_accessibility_divergence.pdf")


pairwise.wilcox.test(abs(ALL$P_est.mean),ALL$class)

### magnitude of divergence SIG DIFF ####
Q <- ALL[ALL$P_qvalue < 0.05 | ALL$H_qvalue < 0.05,] %>%
ggplot(aes(x=class, y=abs(P_est.mean), fill=class)) +
geom_boxplot(notch=TRUE) +
ylim(0,1)

### frequency of sig different
A = nrow(ALL[ALL$class == "inter" & ALL$P_qvalue < 0.05,])/(nrow(ALL[ALL$class == "inter",]))
B = nrow(ALL[ALL$class == "intra" & ALL$P_qvalue < 0.05,])/(nrow(ALL[ALL$class == "intra",]))
C = nrow(ALL[ALL$class == "start" & ALL$P_qvalue < 0.05,])/(nrow(ALL[ALL$class == "start",]))
D = nrow(ALL[ALL$class == "end" & ALL$P_qvalue < 0.05,])/(nrow(ALL[ALL$class == "end",]))

freq_df <- data.frame(
  class = c("inter", "intra", "start", "end"),
  value = c(A, B, C, D)
)

freq_df$class <- factor(freq_df$class,levels = c("start", "end", "inter", "intra"))


Z <- ggplot(freq_df, aes(x=class, y=value)) +
geom_col(aes(fill=class), width=0.2) +
coord_flip() +
ylim(0,1) +
scale_fill_discrete(guide=FALSE) +
ylab("Signficant:non-signficant (q < 0.05)") +
xlab("") +
ggtitle("Proportion of differentially accessible regions") +
theme_main() +
scale_x_discrete(labels=c("txStart","txEnd", "intergenic", "intragenic"))
ggsave(Z, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/ZHR_Z30_Full_results_output_ALL_classes_diff_acess_proportion.pdf")


##### How does variation correlate with accessibility divergence at different regions? ########
## compute number of SNPs within
bedtools -wa -a /Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/ZHR_Z30_Full_results_output_ALL_classes_shortened.bed -b /Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/zhr_z30_SNPs_dm6_genome.bed > temp.bed

df <- read.table("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/temp.bed", header = F)
df_count <- ddply(df,.(V4),nrow)
colnames(df_count) <- c("Paste_locus", "variant_count")
ALL <- join_all(list(ALL, df_count), by = "Paste_locus", type = "full")
ALL$region_length <- ALL$end - ALL$start
ALL$variant_count_length_norm <- (ALL$variant_count/ALL$region_length)*1000

ALL[ALL$class == "intra" & (ALL$P_qvalue < 0.05 | ALL$H_qvalue < 0.05),] %>%
ggplot(aes(x=perc_cis, y=variant_count_length_norm)) +
geom_point()


ALL[ALL$P_qvalue < 0.05 | ALL$H_qvalue < 0.05,] %>%
ggplot(aes(x=variant_count_length_norm,fill=class)) +
geom_density(alpha=0.5)

## calc intragenic regions location
intragenic <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/intragenic_CDS_overlap.bed", header = F)
intragenic$Paste_locus <- paste(intragenic$V1, intragenic$V2, intragenic$V3, sep="_")
intragenic <- intragenic[,c(5,4)]
colnames(intragenic)[2] <- "CDS_overlap"
intragenic$Class <- "NA"
for (i in 1:nrow(intragenic)) {

if (intragenic$CDS_overlap[i] == 0){

	intragenic$Class[i] <- "Intron"

} else if (intragenic$CDS_overlap[i] > 0 & intragenic$CDS_overlap[i] < 1000){

	intragenic$Class[i] <- "Intron_Exon_hyb"

} else if (intragenic$CDS_overlap[i] == 1000){

	intragenic$Class[i] <- "Exon"

}
}
intragenic$CDS_overlap <- NULL

ALL2 <- left_join(ALL, intragenic, by = "Paste_locus") %>% na.omit() %>% unique()

nrow(ALL2[ALL2$Class == "Intron",])
nrow(ALL2[ALL2$Class == "Exon",])
nrow(ALL2[ALL2$Class == "Intron_Exon_hyb",])
