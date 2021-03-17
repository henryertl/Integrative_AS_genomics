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
#######################################################################
###### Generate base plots - % cis accumulation and differences #######
#######################################################################

# Read in Bayes output file
ALL <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_TSIM_ATAC_20min_downsamp_overlap.txt", header = T) %>% na.omit() %>% unique()

## factor sort to have classes in same order in every plot
ALL$class <- factor(ALL$class,levels = c("start", "end", "inter", "intra"))

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
ggsave(R, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/ZHR_TSIM_Full_results_output_ALL_classes_perc_cis.pdf")

pairwise.wilcox.test(ALL$perc_cis,ALL$class)

### magnitude of divergence ALL ####
S <- ggplot(ALL,aes(x=class, y=abs(P_est.mean), fill=class)) +
geom_boxplot(notch=TRUE) +
theme_main() +
ylim(0,1) +
ylab("Estimated accessibility divergence") +
xlab("") +
scale_fill_discrete(guide=FALSE) +
scale_x_discrete(labels=c("txStart","txEnd", "intergenic", "intragenic"))
ggsave(S, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/ZHR_TSIM_Full_results_output_ALL_classes_accessibility_divergence.pdf")


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
ggsave(Z, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/ZHR_TSIM_Full_results_output_ALL_classes_diff_acess_proportion.pdf")


### reinforcing vs opposing
#### only focus on cis-trans changes
##Run classifier for opposing and same
Full_results_output <- ALL

Full_results_output$Direction <- "NA"

for (i in 1:nrow(Full_results_output)) {

if (Full_results_output$Regulatory_class[i] == "Cis_+_Trans,opposing"){

	Full_results_output$Direction[i] <- "Opposing"

} else if (Full_results_output$Regulatory_class[i] == "Cis_*_Trans"){

	Full_results_output$Direction[i] <- "Opposing"

} else if (Full_results_output$Regulatory_class[i] == "Compensatory"){

	Full_results_output$Direction[i] <- "Opposing"

} else if (Full_results_output$Regulatory_class[i] == "Cis_+_Trans,same"){

	Full_results_output$Direction[i] <- "Reinforcing"
}
}

ALL_cis_trans <- Full_results_output


#### create df with reinforcing and opposing numbers for each class
M <- nrow(ALL_cis_trans[ALL_cis_trans$class == "inter" & ALL_cis_trans$Direction == "Opposing",])/nrow(ALL_cis_trans[ALL_cis_trans$class == "inter" & (ALL_cis_trans$Direction == "Opposing" | ALL_cis_trans$Direction == "Reinforcing"),])
N <- nrow(ALL_cis_trans[ALL_cis_trans$class == "intra" & ALL_cis_trans$Direction == "Opposing",])/nrow(ALL_cis_trans[ALL_cis_trans$class == "intra" & (ALL_cis_trans$Direction == "Opposing" | ALL_cis_trans$Direction == "Reinforcing"),])
O <- nrow(ALL_cis_trans[ALL_cis_trans$class == "start" & ALL_cis_trans$Direction == "Opposing",])/nrow(ALL_cis_trans[ALL_cis_trans$class == "start" & (ALL_cis_trans$Direction == "Opposing" | ALL_cis_trans$Direction == "Reinforcing"),])
P <- nrow(ALL_cis_trans[ALL_cis_trans$class == "end" & ALL_cis_trans$Direction == "Opposing",])/nrow(ALL_cis_trans[ALL_cis_trans$class == "end" & (ALL_cis_trans$Direction == "Opposing" | ALL_cis_trans$Direction == "Reinforcing"),])
Q <- nrow(ALL_cis_trans[ALL_cis_trans$class == "inter" & ALL_cis_trans$Direction == "Reinforcing",])/nrow(ALL_cis_trans[ALL_cis_trans$class == "inter" & (ALL_cis_trans$Direction == "Opposing" | ALL_cis_trans$Direction == "Reinforcing"),])
R <- nrow(ALL_cis_trans[ALL_cis_trans$class == "intra" & ALL_cis_trans$Direction == "Reinforcing",])/nrow(ALL_cis_trans[ALL_cis_trans$class == "intra" & (ALL_cis_trans$Direction == "Opposing" | ALL_cis_trans$Direction == "Reinforcing"),])
S <- nrow(ALL_cis_trans[ALL_cis_trans$class == "start" & ALL_cis_trans$Direction == "Reinforcing",])/nrow(ALL_cis_trans[ALL_cis_trans$class == "start" & (ALL_cis_trans$Direction == "Opposing" | ALL_cis_trans$Direction == "Reinforcing"),])
U <- nrow(ALL_cis_trans[ALL_cis_trans$class == "end" & ALL_cis_trans$Direction == "Reinforcing",])/nrow(ALL_cis_trans[ALL_cis_trans$class == "end" & (ALL_cis_trans$Direction == "Opposing" | ALL_cis_trans$Direction == "Reinforcing"),])

freq_df <- data.frame(
  class = c("inter", "intra", "start", "end"),
  opposing = c(M, N, O, P),
  reinforcing = c(Q, R, S, U)
)

L <- melt(freq_df) %>%
ggplot(aes(x=class, y=value)) +
geom_col(aes(fill=variable), width = 0.2) +
coord_flip() +
ylim(0,1) +
ylab("opposing:reinforcing cis-trans changes") +
xlab("") +
geom_hline(yintercept=0.5, linetype="dashed", color = "grey") +
theme_main() +
scale_x_discrete(labels=c("txStart","txEnd", "intergenic", "intragenic"))
ggsave(L, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/ZHR_TSIM_Full_results_output_ALL_classes_oppos_reinforc_proportion.pdf")


####### How does variation correlate with accessibility divergence at different regions? ########
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
