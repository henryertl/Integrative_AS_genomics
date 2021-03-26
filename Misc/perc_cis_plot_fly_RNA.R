# load libraries and functions
library(stringr)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(rstatix)

setwd("/Users/wittkopp_member/Code")

theme_main <- function() {
  theme_bw() +
  theme(
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  axis.text = element_text(size = 30),
  axis.title = element_text(size = 40),
  strip.text = element_text(size = 30),
  legend.text= element_text(size = 30),
  legend.title = element_text(size = 20),
  plot.title = element_text(size = 25, face = "bold")

)
}

#######################################################################
###### Generate base plots - % cis accumulation and differences #######
#######################################################################

# Read in Bayes output files
TSIM <- read.delim("Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_TSIM_RNA_20min.txt", header = T)
Z30 <- read.delim("Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_1000max.txt", header = T)

# Merge from both comparisons
Full_results_comb <- merge(Z30, TSIM, by = "Paste_locus")
Full_results_comb <- Full_results_comb %>%
  mutate(Gene = str_sub(Paste_locus, 7, -1))
Full_results_comb <- Full_results_comb[-c(1)]
Full_results_comb$perc_cis_diff <- log2(Full_results_comb$perc_cis.y/Full_results_comb$perc_cis.x)

write.table(Full_results_comb, file="Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/RNA_seq/ZHR_Z30_TSIM_Full_results_output.txt", sep = "\t", row.names = F, quote = F)

Full_results_comb <- read.delim("/AS_ATAC_RNA_2020_10_1/RNA_seq/ZHR_Z30_TSIM_Full_results_output.txt", header = T)

# Plot % cis across species
Full_results_comb_div <- Full_results_comb[Full_results_comb$P_qvalue.x < 0.05 | Full_results_comb$H_qvalue.x < 0.05 | Full_results_comb$P_qvalue.y < 0.05 | Full_results_comb$H_qvalue.y < 0.05,]
Perc_cis <- cbind(Full_results_comb_div $perc_cis.x, Full_results_comb_div $perc_cis.y)
colnames(Perc_cis) <- c("Mel_Mel", "Mel_Sim")
Perc_cis_melt <- melt(Perc_cis)

R  <- ggplot(Perc_cis_melt, aes(x=Var2, y=value, fill=Var2)) +
geom_boxplot(notch=TRUE) +
scale_x_discrete(labels=c("Mel-Mel","Mel-Sim")) +
xlab("") +
ylab("% cis") +
scale_fill_discrete(guide=FALSE) +
theme_main() +
stat_compare_means(method = "wilcox.test")
  ggsave(R, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Figures/perc_cis_within_and_between_RNA.pdf", width = 10, height = 15)

# Plot distribution of ∆ cis by divergence
Perc_cis <- as.data.frame(Perc_cis)
Perc_cis$perc_cis_diff <- log2(Perc_cis$Mel_Sim/Perc_cis$Mel_Mel)

S <- ggplot(Perc_cis, aes(x=perc_cis_diff)) +
geom_density(fill = "#E69F00") +
theme_main() +
xlab("Change in cis by divergence") +
geom_vline(xintercept=0)
  ggsave(S, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Figures/perc_cis_divergence_distribution_within_and_between_RNA.pdf", width = 10, height = 15)

#######################################################################
###### compare % cis accumulation for different characteristics #######
#######################################################################

########################### X VS AUTOSOMES ############################

# Append relevant data to df
# Get gene coordinates (need to rip chromosomes from other file)
df_exon_counts <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/ZHR_TSIM_genic_counts_dm6.bed", header = F)
gene_coords <- df_exon_counts[, c("V1", "V4")]
colnames(gene_coords) <- c("chrom", "Gene")
Full_results_comb_coords <- join_all(list(gene_coords, y=Full_results_comb), by = 'Gene', type = 'right') %>% unique()

Full_results_comb_coords <- Full_results_comb_coords %>% drop_na("chrom")

# Assign X or Autosome category
Full_results_comb_coords$Auto_X <- "NA"

for (i in 1:nrow(Full_results_comb_coords)) {

if (Full_results_comb_coords$chrom[i] == "chrX"){

	Full_results_comb_coords$Auto_X[i] <- "X"

} else if (Full_results_comb_coords$chrom[i] == "chr2L"){

	Full_results_comb_coords$Auto_X[i] <- "Autosome"

} else if (Full_results_comb_coords$chrom[i] == "chr2R"){

	Full_results_comb_coords$Auto_X[i] <- "Autosome"

} else if (Full_results_comb_coords$chrom[i] == "chr3L"){

	Full_results_comb_coords$Auto_X[i] <- "Autosome"

} else if (Full_results_comb_coords$chrom[i] == "chr3R"){

	Full_results_comb_coords$Auto_X[i] <- "Autosome"

} else if (Full_results_comb_coords$chrom[i] == "chr4"){

Full_results_comb_coords$Auto_X[i] <- "Autosome"

}
}

# IS THERE A FASTER X EFFECT IN THIS DATA?

## Plot distribution of estimated expression differences (P_est.mean)
ggplot(Full_results_comb_coords, aes(x=Auto_X, y=abs(P_est.mean.y), fill = Auto_X))+
geom_violin() +
geom_boxplot(notch=TRUE, width=0.2) +
#ylim(-2.5,2.5) +
xlab("") +
ylab("Estimated expression difference") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
stat_compare_means(method="wilcox.test")

## Plot ratio of genes with parental sig difference to non-sig difference
A <- nrow(Full_results_comb_coords[Full_results_comb_coords$P_p_value.y < 0.05 & Full_results_comb_coords$Auto_X == "X",])
B <- nrow(Full_results_comb_coords[Full_results_comb_coords$P_p_value.y >= 0.05 & Full_results_comb_coords$Auto_X == "X",])
C <- A+B
D <- nrow(Full_results_comb_coords[Full_results_comb_coords$P_p_value.y < 0.05 & Full_results_comb_coords$Auto_X == "Autosome",])
E <- nrow(Full_results_comb_coords[Full_results_comb_coords$P_p_value.y >= 0.05 & Full_results_comb_coords$Auto_X == "Autosome",])
G <- D + E

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("X", 2), rep("Autosome", 2))
twobytwo[,2] <- c(A/C,B/C,D/G,E/G)
twobytwo[,3] <- c("Sig","Non_sig","Sig","Non_sig")
colnames(twobytwo) <- c("Comparison", "Proportion", "Chromosome_class")

ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Chromosome_class)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  ggtitle("Proportion of genes with significant differential expression")

# FET = 0.0055

## Plot spearman correlation between replicates
x <- Full_results_comb_coords[Full_results_comb_coords$Auto_X == "X",]
auto <- Full_results_comb_coords[Full_results_comb_coords$Auto_X == "Autosome",]

P1_cor_X <- x[, c("P1_1.y", "P2_1.y")]
P1_cor_Auto <- auto[, c("P1_1.y", "P2_1.y")]
X1 <- cor(P1_cor_X, method = "spearman")
A1 <- cor(P1_cor_Auto, method = "spearman")

P2_cor_X <- x[, c("P1_2.y", "P2_2.y")]
P2_cor_Auto <- auto[, c("P1_2.y", "P2_2.y")]
X2 <- cor(P2_cor_X, method = "spearman")
A2 <- cor(P2_cor_Auto, method = "spearman")

P3_cor_X <- x[, c("P1_3.y", "P2_3.y")]
P3_cor_Auto <- auto[, c("P1_3.y", "P2_3.y")]
X3 <- cor(P3_cor_X, method = "spearman")
A3 <- cor(P3_cor_Auto, method = "spearman")

P4_cor_X <- x[, c("P1_4.y", "P2_4.y")]
P4_cor_Auto <- auto[, c("P1_4.y", "P2_4.y")]
X4 <- cor(P4_cor_X, method = "spearman")
A4 <- cor(P4_cor_Auto, method = "spearman")

P5_cor_X <- x[, c("P1_5.y", "P2_5.y")]
P5_cor_Auto <- auto[, c("P1_5.y", "P2_5.y")]
X5 <- cor(P5_cor_X, method = "spearman")
A5 <- cor(P5_cor_Auto, method = "spearman")

P6_cor_X <- x[, c("P1_6.y", "P2_6.y")]
P6_cor_Auto <- auto[, c("P1_6.y", "P2_6.y")]
X6 <- cor(P6_cor_X, method = "spearman")
A6 <- cor(P6_cor_Auto, method = "spearman")


X_A_spear_corr <- NULL
X_A_spear_corr$X <- c(X1[1,2], X2[1,2], X3[1,2], X4[1,2], X5[1,2], X6[1,2])
X_A_spear_corr$A <- c(A1[1,2], A2[1,2], A3[1,2], A4[1,2], A5[1,2], A6[1,2])
X_A_spear_corr <- as.data.frame(X_A_spear_corr)
X_A_spear_corr$X <- 1 - X_A_spear_corr$X
X_A_spear_corr$A <- 1 - X_A_spear_corr$A
melt(X_A_spear_corr) %>%
ggpaired(x='variable', y='value') +
stat_compare_means(paired = TRUE) +
ylab("1 - rho (Spearmans corr coefficient)") +
xlab("")

# IS PERCENT CIS CHANGE GREATER FOR X THAN FOR AUTOSOMES?
# Plot ∆percent cis difference: X versus Autosomes
Full_results_comb_coords[Full_results_comb_coords$P_p_value.x < 0.05 | Full_results_comb_coords$H_p_value.x < 0.05,] %>%
ggplot(aes(x=Auto_X, y=perc_cis_diff, fill = Auto_X))+
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
xlab("") +
ylab("Percent cis difference") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
stat_compare_means(method="wilcox.test")
  ggsave(B, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Figures/Expression_difference_X_Auto.pdf", width = 5, height = 15)


######################### EXPRESSION DIVERGENCE ##########################

# Filter only genes with expression divergence in either parents or hybrid
Full_results_comb_exp_div <- Full_results_comb[Full_results_comb$P_p_value.y < 0.05 | Full_results_comb$H_p_value.y < 0.05,] %>% unique() %>% drop_na("P_est.mean.y") %>% drop_na("P1_1.y")

# Plot distribution of expression divergence just to get a sense of the data
ggplot(Full_results_comb_exp_div, aes(abs(P_est.mean.y))) +
geom_density(fill = "#E69F00") +
theme_main() +
xlab("Estimated expression divergence Mel-Sim")

## Get summary stats on expression for P1_1.x
stat_sum <- summary(abs(Full_results_comb_exp_div$P_est.mean.y))

## Assign expression change quartile to each gene
Full_results_comb_exp_div$Expression_change_quantile <- "NA"

for (i in 1:nrow(Full_results_comb_exp_div)) {

if (Full_results_comb_exp_div$P_est.mean.y[i] <= as.numeric(stat_sum[2])){

	Full_results_comb_exp_div$Expression_change_quantile[i] <- "1Q"

} else if (Full_results_comb_exp_div$P_est.mean.y[i] >= as.numeric(stat_sum[2]) & Full_results_comb_exp_div$P_est.mean.y[i] <= as.numeric(stat_sum[3])){

	Full_results_comb_exp_div$Expression_change_quantile[i] <- "2Q"

} else if (Full_results_comb_exp_div$P_est.mean.y[i] >= as.numeric(stat_sum[3]) & Full_results_comb_exp_div$P_est.mean.y[i] <= as.numeric(stat_sum[5])){

	Full_results_comb_exp_div$Expression_change_quantile[i] <- "3Q"

} else if (Full_results_comb_exp_div$P_est.mean.y[i] >= as.numeric(stat_sum[5]) & Full_results_comb_exp_div$P_est.mean.y[i] <= as.numeric(stat_sum[6])){

	Full_results_comb_exp_div$Expression_change_quantile[i] <- "4Q"

}
}

# Assign top 10% expression change or other
dec_stat <- quantile(abs(Full_results_comb_exp_div$P_est.mean.y), prob = seq(0,1, length = 11), type = 5)

Full_results_comb_exp_div$Expression_change_dec <- "NA"

for (i in 1:nrow(Full_results_comb_exp_div)) {

if (abs(Full_results_comb_exp_div$P_est.mean.y[i]) < as.numeric(dec_stat[10])){

	Full_results_comb_exp_div$Expression_change_dec[i] <- "0_89"

} else if (abs(Full_results_comb_exp_div$P_est.mean.y[i]) >= as.numeric(dec_stat[10])){

	Full_results_comb_exp_div$Expression_change_dec[i] <- "90_100"

}
}

######
# Filter only genes with expression divergence in either parents or hybrid
Full_results_comb_exp_div <- Full_results_comb[Full_results_comb$P_p_value.x < 0.05 | Full_results_comb$H_p_value.x < 0.05,] %>% unique() %>% drop_na("P_est.mean.x") %>% drop_na("P1_1.x")


## Get summary stats on expression for P1_1.x
stat_sum <- summary(abs(Full_results_comb_exp_div$P_est.mean.x))

## Assign expression change quartile to each gene
Full_results_comb_exp_div$Expression_change_quantile <- "NA"

for (i in 1:nrow(Full_results_comb_exp_div)) {

if (Full_results_comb_exp_div$P_est.mean.x[i] <= as.numeric(stat_sum[2])){

	Full_results_comb_exp_div$Expression_change_quantile[i] <- "1Q"

} else if (Full_results_comb_exp_div$P_est.mean.x[i] >= as.numeric(stat_sum[2]) & Full_results_comb_exp_div$P_est.mean.x[i] <= as.numeric(stat_sum[3])){

	Full_results_comb_exp_div$Expression_change_quantile[i] <- "2Q"

} else if (Full_results_comb_exp_div$P_est.mean.x[i] >= as.numeric(stat_sum[3]) & Full_results_comb_exp_div$P_est.mean.x[i] <= as.numeric(stat_sum[5])){

	Full_results_comb_exp_div$Expression_change_quantile[i] <- "3Q"

} else if (Full_results_comb_exp_div$P_est.mean.x[i] >= as.numeric(stat_sum[5]) & Full_results_comb_exp_div$P_est.mean.x[i] <= as.numeric(stat_sum[6])){

	Full_results_comb_exp_div$Expression_change_quantile[i] <- "4Q"

}
}

# Assign top 10% expression change or other
dec_stat <- quantile(abs(Full_results_comb_exp_div$P_est.mean.x), prob = seq(0,1, length = 11), type = 5)

Full_results_comb_exp_div$Expression_change_dec <- "NA"

for (i in 1:nrow(Full_results_comb_exp_div)) {

if (abs(Full_results_comb_exp_div$P_est.mean.x[i]) < as.numeric(dec_stat[10])){

	Full_results_comb_exp_div$Expression_change_dec[i] <- "0_89"

} else if (abs(Full_results_comb_exp_div$P_est.mean.x[i]) >= as.numeric(dec_stat[10])){

	Full_results_comb_exp_div$Expression_change_dec[i] <- "90_100"

}
}
######

# Plot
## ARE TOP 10% OF DIFFERENTIALLY EXPRESSED GENES MORE % CIS DIFF THAN THE REST?
A <- ggplot(Full_results_comb_exp_div, aes(x = Expression_change_dec, y = perc_cis_diff, fill = Expression_change_dec))+
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
#ylim(-2.5,2.5) +
xlab("Expression divergence quartile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
stat_compare_means(method="wilcox.test")

## RELATIONSHIP BETWEEN EXPRESSION DIVERGENCE AND PERCENT CIS DIFFERENCE
B <-ggplot(Full_results_comb_exp_div, aes(x = Expression_change_quantile, y = perc_cis_diff, fill = Expression_change_quantile))+
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.3) +
#ylim(-5,5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE)

ggplot(Full_results_comb_exp_div, aes(y=abs(P_est.mean.x), x=perc_cis_diff)) +
geom_point() +
xlim(-2,2) +
ylim(0,2)

ggplot(Full_results_comb) +
geom_violin(aes(x=abs(P_est.mean.x), x=abs(P_est.mean.y)))


summary(abs(Full_results_comb$P_est.mean.x))
summary(abs(Full_results_comb$P_est.mean.y))



## DO THE SAME WITHOUT X CHROMOSOME GENES

Full_results_comb_exp_div_coords <- left_join(Full_results_comb_exp_div, Full_results_comb_coords, by = "Gene")

A <- Full_results_comb_exp_div_coords[Full_results_comb_exp_div_coords$Auto_X == "Autosome",] %>%
ggplot(aes(x = Expression_change_quantile, y = perc_cis_diff.y, fill = Expression_change_quantile))+
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.3) +
#ylim(-5,5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE)

Full_results_comb_exp_div_coords[Full_results_comb_exp_div_coords$Auto_X == "Autosome",] %>%
ggplot(aes(x = Expression_change_dec.y, y = perc_cis_diff.y, fill = Expression_change_dec.y))+
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
#ylim(-2.5,2.5) +
xlab("Expression divergence quartile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
stat_compare_means(method="wilcox.test")


################# EXPRESSION LEVEL ####################

## Get summary stats on expression for P1_1.x
stat_sum <- summary(abs(Full_results_comb_exp_div$P1_1.x))

ggplot(Full_results_comb_exp_div, aes(x=P1_1.x))+
geom_density()

## Assign expression change quartile to each gene
Full_results_comb_exp_div$Expression_level_quantile <- "NA"

for (i in 1:nrow(Full_results_comb_exp_div)) {

if (Full_results_comb_exp_div$P1_1.x[i] <= as.numeric(stat_sum[2])){

	Full_results_comb_exp_div$Expression_level_quantile[i] <- "1Q"

} else if (Full_results_comb_exp_div$P1_1.x[i] >= as.numeric(stat_sum[2]) & Full_results_comb_exp_div$P1_1.x[i] <= as.numeric(stat_sum[3])){

	Full_results_comb_exp_div$Expression_level_quantile[i] <- "2Q"

} else if (Full_results_comb_exp_div$P1_1.x[i] >= as.numeric(stat_sum[3]) & Full_results_comb_exp_div$P1_1.x[i] <= as.numeric(stat_sum[5])){

	Full_results_comb_exp_div$Expression_level_quantile[i] <- "3Q"

} else if (Full_results_comb_exp_div$P1_1.x[i] >= as.numeric(stat_sum[5]) & Full_results_comb_exp_div$P1_1.x[i] <= as.numeric(stat_sum[6])){

	Full_results_comb_exp_div$Expression_level_quantile[i] <- "4Q"

}
}

# Assign top 10% expression change or other
dec_stat <- quantile(abs(Full_results_comb_exp_div$P1_1.y), prob = seq(0,1, length = 11), type = 5)

Full_results_comb_exp_div$Expression_level_quantile_dec <- "NA"

for (i in 1:nrow(Full_results_comb_exp_div)) {

if (abs(Full_results_comb_exp_div$P1_1.y[i]) < as.numeric(dec_stat[10])){

	Full_results_comb_exp_div$Expression_level_quantile_dec[i] <- "0_89"

} else if (abs(Full_results_comb_exp_div$P1_1.y[i]) >= as.numeric(dec_stat[10])){

	Full_results_comb_exp_div$Expression_level_quantile_dec[i] <- "90_100"

}
}


# Plot
## ARE TOP 10% OF EXPRESSED GENES MORE % CIS DIFF THAN THE REST?
R <- ggplot(Full_results_comb_exp_div, aes(x = Expression_level_quantile_dec, y = perc_cis_diff, fill = Expression_level_quantile_dec))+
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
#ylim(-2.5,2.5) +
xlab("Expression divergence quartile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
stat_compare_means(method="wilcox.test")

## RELATIONSHIP BETWEEN EXPRESSION LEVEL AND PERCENT CIS DIFFERENCE
S <- ggplot(Full_results_comb_exp_div, aes(x = Expression_level_quantile, y = perc_cis_diff, fill = Expression_level_quantile))+
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
#ylim(-5,5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE)

S

ggplot(Full_results_comb_exp_div) +
geom_point(aes(x=perc_cis.y, y=abs(P_est.mean.y))


######################### NETWORK CHARACTERISTICS ##########################

# Network properties
## Read in network info
network_props <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/ZHR.SIM.INDEX.NEW.txt", header = T)
network_props <- network_props[,c(1,2,4)]

ID_conversion_table <- read.delim("/Users/henryertl/Documents/Bing_network_ID_conversion.txt", header = T)
colnames(ID_conversion_table) <- c("geneID", "Gene")

network_props <- full_join(network_props, ID_conversion_table, by = "geneID")

Full_results_comb_coords_network <- left_join(Full_results_comb, network_props, by = "Gene") %>% as.data.frame()

## Indegree
## Assign inDEGREE quartile to each gene
Full_results_comb_coords_network_indegree <- Full_results_comb_coords_network %>% drop_na("inDEGREE")

stat_sum <- summary(abs(Full_results_comb_coords_network_indegree$inDEGREE))

Full_results_comb_coords_network_indegree$in_deg_quantile <- "NA"

for (i in 1:nrow(Full_results_comb_coords_network_indegree)) {

if (Full_results_comb_coords_network_indegree$inDEGREE[i] <= as.numeric(stat_sum[2])){

	Full_results_comb_coords_network_indegree$in_deg_quantile[i] <- "1Q"

} else if (Full_results_comb_coords_network_indegree$inDEGREE[i] >= as.numeric(stat_sum[2]) & Full_results_comb_coords_network_indegree$inDEGREE[i] <= as.numeric(stat_sum[3])){

	Full_results_comb_coords_network_indegree$in_deg_quantile[i] <- "2Q"

} else if (Full_results_comb_coords_network_indegree$inDEGREE[i] >= as.numeric(stat_sum[3]) & Full_results_comb_coords_network_indegree$inDEGREE[i] <= as.numeric(stat_sum[5])){

	Full_results_comb_coords_network_indegree$in_deg_quantile[i] <- "3Q"

} else if (Full_results_comb_coords_network_indegree$inDEGREE[i] >= as.numeric(stat_sum[5]) & Full_results_comb_coords_network_indegree$inDEGREE[i] <= as.numeric(stat_sum[6])){

	Full_results_comb_coords_network_indegree$in_deg_quantile[i] <- "4Q"

}
}

dec_stat <- quantile(Full_results_comb_coords_network_indegree$inDEGREE, prob = seq(0,1, length = 11), type = 5)

Full_results_comb_coords_network_indegree$in_deg_dec <- "NA"

for (i in 1:nrow(Full_results_comb_coords_network_indegree)) {

if (Full_results_comb_coords_network_indegree$inDEGREE[i] < as.numeric(dec_stat[10])){

	Full_results_comb_coords_network_indegree$in_deg_dec[i] <- "0_89"

} else if (Full_results_comb_coords_network_indegree$inDEGREE[i] >= as.numeric(dec_stat[10])){

	Full_results_comb_coords_network_indegree$in_deg_dec[i] <- "90_100"

}
}



L <- ggplot(Full_results_comb_coords_network_indegree, aes(x = in_deg_quantile, y = perc_cis_diff, fill = in_deg_quantile))+
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
#ylim(-2.5,2.5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE)

Z <- ggplot(Full_results_comb_coords_network_indegree, aes(x = in_deg_dec, y = perc_cis_diff, fill = in_deg_dec))+
#geom_violin() +
#geom_violin() +
geom_boxplot(notch=TRUE) +
ylim(-5,5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
stat_compare_means(method="wilcox.test")




## Outdegree
## Assign outDEGREE quartile to each gene
Full_results_comb_coords_network_outdegree <- Full_results_comb_coords_network[Full_results_comb_coords_network$outDEGREE != 0,] %>% drop_na("outDEGREE")

stat_sum <- summary(Full_results_comb_coords_network_outdegree$outDEGREE)

Full_results_comb_coords_network_outdegree$in_deg_quantile <- "NA"

for (i in 1:nrow(Full_results_comb_coords_network_outdegree)) {

if (Full_results_comb_coords_network_outdegree$outDEGREE[i] <= as.numeric(stat_sum[2])){

	Full_results_comb_coords_network_outdegree$out_deg_quantile[i] <- "1Q"

} else if (Full_results_comb_coords_network_outdegree$outDEGREE[i] >= as.numeric(stat_sum[2]) & Full_results_comb_coords_network_outdegree$outDEGREE[i] <= as.numeric(stat_sum[3])){

	Full_results_comb_coords_network_outdegree$out_deg_quantile[i] <- "2Q"

} else if (Full_results_comb_coords_network_outdegree$outDEGREE[i] >= as.numeric(stat_sum[3]) & Full_results_comb_coords_network_outdegree$outDEGREE[i] <= as.numeric(stat_sum[5])){

	Full_results_comb_coords_network_outdegree$out_deg_quantile[i] <- "3Q"

} else if (Full_results_comb_coords_network_outdegree$outDEGREE[i] >= as.numeric(stat_sum[5]) & Full_results_comb_coords_network_outdegree$outDEGREE[i] <= as.numeric(stat_sum[6])){

	Full_results_comb_coords_network_outdegree$out_deg_quantile[i] <- "4Q"

}
}

dec_stat <- quantile(Full_results_comb_coords_network_outdegree$outDEGREE, prob = seq(0,1, length = 11), type = 5)

Full_results_comb_coords_network_outdegree$out_deg_dec <- "NA"

for (i in 1:nrow(Full_results_comb_coords_network_outdegree)) {

if (Full_results_comb_coords_network_outdegree$outDEGREE[i] < 102){

	Full_results_comb_coords_network_outdegree$out_deg_dec[i] <- "0_89"

} else if (Full_results_comb_coords_network_outdegree$outDEGREE[i] >= 102){

	Full_results_comb_coords_network_outdegree$out_deg_dec[i] <- "90_100"

}
}



ggplot(Full_results_comb_coords_network_outdegree, aes(x = out_deg_quantile, y = perc_cis_diff, fill = out_deg_quantile))+
#geom_violin() +
geom_boxplot(notch=FALSE, width = 0.2) +
#ylim(-2.5,2.5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE)

ggplot(Full_results_comb_coords_network_outdegree, aes(x = out_deg_dec, y = perc_cis_diff, fill = out_deg_dec))+
geom_boxplot() +
ylim(-5,5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
stat_compare_means(method="wilcox.test")

# outdegree + indegree

Full_results_comb_coords_network_indegree_outdegree <- Full_results_comb_coords_network %>% drop_na("inDEGREE") %>% drop_na("outDEGREE") %>% unique()

# Get distributions of both in and out degree

ggplot(Full_results_comb_coords_network_indegree_outdegree, aes(inDEGREE)) +
geom_density()

ggplot(Full_results_comb_coords_network_indegree_outdegree, aes(outDEGREE)) +
geom_density() +
xlim(0,50)

# Classify as highly connected or not

Full_results_comb_coords_network_indegree_outdegree$connectivity <- "NO"

for (i in 1:nrow(Full_results_comb_coords_network_indegree_outdegree)) {

if (Full_results_comb_coords_network_indegree_outdegree$outDEGREE[i] > 20 & Full_results_comb_coords_network_indegree_outdegree$outDEGREE[i] > 20){

	Full_results_comb_coords_network_indegree_outdegree$connectivity[i] <- "YES"

}
}

B <- ggplot(Full_results_comb_coords_network_indegree_outdegree, aes(x = connectivity, y = perc_cis_diff, fill = connectivity))+
geom_boxplot() +
#ylim(-5,5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
stat_compare_means(method="wilcox.test")

B
####### Extra Analyses ####


Full_results_comb_coords[Full_results_comb_coords$outDEGREE != 0,] %>%

ggplot(Full_results_comb_coords, aes(inDEGREE)) +
geom_density()

+
xlim(0,5)


Full_results_comb_coords$in_deg_bin <- "NA"

for (i in 1:nrow(Full_results_comb_coords)) {

if (Full_results_comb_coords$inDEGREE[i] <= 25){

	Full_results_comb_coords$in_deg_bin[i] <- "1"

} else if (Full_results_comb_coords$inDEGREE[i] >= 25 & Full_results_comb_coords$inDEGREE[i] <= 50){

	Full_results_comb_coords$in_deg_bin[i] <- "2"

} else if (Full_results_comb_coords$inDEGREE[i] >= 50 & Full_results_comb_coords$inDEGREE[i] <= 75){

	Full_results_comb_coords$in_deg_bin[i] <- "3"

} else if (Full_results_comb_coords$inDEGREE[i] >= 75 & Full_results_comb_coords$inDEGREE[i] <= 100){

	Full_results_comb_coords$in_deg_bin[i] <- "4"

} else if (Full_results_comb_coords$inDEGREE[i] >= 100 & Full_results_comb_coords$inDEGREE[i] <= 125){

	Full_results_comb_coords$in_deg_bin[i] <- "5"

} else if (Full_results_comb_coords$inDEGREE[i] >= 125 & Full_results_comb_coords$inDEGREE[i] <= 150){

  Full_results_comb_coords$in_deg_bin[i] <- "6"

} else if (Full_results_comb_coords$inDEGREE[i] >= 150 & Full_results_comb_coords$inDEGREE[i] <= 175){

  Full_results_comb_coords$in_deg_bin[i] <- "7"
}
}

stat.test <- Full_results_comb_coords %>%
  wilcox_test(perc_cis_diff ~ in_deg_dec) %>%
  add_significance()
stat.test


quantile(Full_results_comb_coords$inDEGREE, prob = seq(0,1, length = 11), type = 5)
ggplot(Full_results_comb_coords, aes(x = in_deg_bin, y = perc_cis_diff, fill = in_deg_bin))+
#geom_violin() +
geom_boxplot() +
#ylim(-5,5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
geom_smooth(method="loess", aes(group=1))


euclidean <- function(a, b) sqrt(sum((a - b)^2))

euclidian(Full_results_comb_coords)



library(rtracklayer)
export.bw(coverage(greenleaf),'greenleaf.bw')
