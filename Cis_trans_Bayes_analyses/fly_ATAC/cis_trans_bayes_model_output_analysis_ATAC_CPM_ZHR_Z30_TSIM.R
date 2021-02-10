# load libraries and functions
library(stringr)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(magrittr)
library(cowplot)

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
TSIM <- read.delim("~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/ZHR_Z30_TSIM_universal_peakset_outputs/Full_results_output_ZHR_TSIM_ATAC_20min.txt", header = T)
Z30 <- read.delim("~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/ZHR_Z30_TSIM_universal_peakset_outputs/Full_results_output_ZHR_Z30_ATAC_20min.txt", header = T)

# Merge from both comparisons
Full_results_comb <- merge(Z30, TSIM, by = "Paste_locus") %>% unique()
write.table(Full_results_comb, file="~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/ZHR_Z30_TSIM_universal_peakset_outputs/Full_results_output_ZHR_Z30_TSIM_ATAC_20min.txt", sep = "\t", row.names = F, quote = F)

# Plot % cis across species
Full_results_comb_div <- Full_results_comb[Full_results_comb$Regulatory_class.x != "Conserved/Ambiguous" & Full_results_comb$Regulatory_class.y != "Conserved/Ambiguous",]
Perc_cis <- cbind(Full_results_comb_div $perc_cis.x, Full_results_comb_div $perc_cis.y)
colnames(Perc_cis) <- c("Mel_Mel", "Mel_Sim")
Perc_cis_melt <- melt(Perc_cis)

R  <- ggplot(Perc_cis_melt, aes(x=Var2, y=value, fill=Var2)) +
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
scale_x_discrete(labels=c("Within species","Between species")) +
xlab("") +
ylab("% cis") +
scale_fill_discrete(guide=FALSE) +
theme_main() +
stat_compare_means(method = "wilcox.test")
  ggsave(R, file = "~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/perc_cis_within_and_between_ATAC.pdf", width = 10, height = 15)

# Define as coding non-coding
## Isolate gene coordinates from ATAC peaks
ATAC_coords <- Full_results_comb[,c(2:4,1)]
colnames(ATAC_coords) <- c("chrom", "start", "end", "Paste_coords")
#write.table(ATAC_coords, file="~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Full_results_output_ZHR_Z30_TSIM_ATAC_20min_peaks.bed", sep = "\t", row.names = F, quote = F)

## Read in gene coordinates
coding_coords <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Full_results_output_ZHR_Z30_TSIM_ATAC_20min_peaks_CDS.bed", header = F)
colnames(coding_coords) <- c("chrom","start", "end", "Paste_locus")
non_cod_dist_coords <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Full_results_output_ZHR_Z30_TSIM_ATAC_20min_peaks_non_coding_NON1kBUP.bed", header = F)
colnames(non_cod_dist_coords) <- c("chrom","start", "end", "Paste_locus")
non_cod_prox_coords <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Full_results_output_ZHR_Z30_TSIM_ATAC_20min_peaks_non_coding_1kBUP.bed", header = F)
colnames(non_cod_prox_coords) <- c("chrom","start", "end", "Paste_locus")

## Join and subset by group
Full_results_comb_coding <- merge(x=Full_results_comb, y=coding_coords, by = "Paste_locus")
Full_results_comb_coding$class <- "Coding"
Full_results_comb_non_coding_dist <- merge(x=Full_results_comb, y=non_cod_dist_coords, by = "Paste_locus")
Full_results_comb_non_coding_dist$class <- "Distal"
Full_results_comb_non_coding_prox <- merge(x=Full_results_comb, y=non_cod_prox_coords, by = "Paste_locus") %>% unique()
Full_results_comb_non_coding_prox$class <- "Proximal"
Full_results_comb_non_coding <- rbind(Full_results_comb_non_coding_dist,Full_results_comb_non_coding_prox)


Full_results_comb <- rbind(Full_results_comb_coding,Full_results_comb_non_coding_dist,Full_results_comb_non_coding_prox)
Full_results_comb$abs_P_est.mean.y <- abs(Full_results_comb$P_est.mean.y)
Full_results_comb$abs_P_est.mean.x <- abs(Full_results_comb$P_est.mean.x)

## PLOT VARIATION ACROSS GROUPS
## Estimated accessibility change
A <- ggplot(Full_results_comb, aes(y=abs_P_est.mean.x, x=class, fill=class)) +
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
theme_main() +
ylab("Estimated CA change") +
xlab("")  +
scale_fill_discrete(guide=FALSE) +
ggtitle("Estimated CA change within species comparison") +
ylim(0,2)
ggsave(A, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Est_CA_change_by_group_within.pdf", width = 10, height = 15)

B <- ggplot(Full_results_comb, aes(y=abs_P_est.mean.y, x=class, fill=class)) +
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
theme_main() +
ylab("Estimated CA change") +
xlab("")  +
scale_fill_discrete(guide=FALSE) +
ggtitle("Estimated CA change between species comparison") +
ylim(0,2)
ggsave(B, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Est_CA_change_by_group_between.pdf", width = 10, height = 15)


stat.test <- Full_results_comb %>%
  wilcox_test(abs_P_est.mean.y ~ class)

stat.test <- Full_results_comb %>%
  wilcox_test(abs_P_est.mean.x ~ class)

## Plot ratio of regions with parental sig difference to non-sig difference
A <- nrow(Full_results_comb[Full_results_comb$P_qvalue.y < 0.05 & Full_results_comb$class == "Coding",])
B <- nrow(Full_results_comb[Full_results_comb$P_qvalue.y >= 0.05 & Full_results_comb$class == "Coding",])
C <- A+B
D <- nrow(Full_results_comb[Full_results_comb$P_qvalue.y < 0.05 & Full_results_comb$class == "Proximal",])
E <- nrow(Full_results_comb[Full_results_comb$P_qvalue.y >= 0.05 & Full_results_comb$class == "Proximal",])
G <- D+E
H <- nrow(Full_results_comb[Full_results_comb$P_qvalue.y < 0.05 & Full_results_comb$class == "Distal",])
I <- nrow(Full_results_comb[Full_results_comb$P_qvalue.y >= 0.05 & Full_results_comb$class == "Distal",])
J <- H+I

threebythree <- matrix(ncol = 3, nrow = 6) %>% as.data.frame()
threebythree[,1] <- c(rep("Coding", 2), rep("Proximal", 2), rep("Distal", 2))
threebythree[,2] <- c(A/C,B/C,D/G,E/G,H/J,I/J)
threebythree[,3] <- c("Sig","Non_sig","Sig","Non_sig","Sig","Non_sig")
colnames(threebythree) <- c("Comparison", "Proportion", "Class")

Z <- ggplot(data = threebythree, aes(x=Comparison, y=Proportion, fill=Class)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  ggtitle("Sig vs non-sig CA - between species ")
  ggsave(Z, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Sig_vs_nonSig_CA_change_by_group_between.pdf", width = 10, height = 15)



A <- nrow(Full_results_comb[Full_results_comb$P_qvalue.x < 0.05 & Full_results_comb$class == "Coding",])
B <- nrow(Full_results_comb[Full_results_comb$P_qvalue.x >= 0.05 & Full_results_comb$class == "Coding",])
C <- A+B
D <- nrow(Full_results_comb[Full_results_comb$P_qvalue.x < 0.05 & Full_results_comb$class == "Proximal",])
E <- nrow(Full_results_comb[Full_results_comb$P_qvalue.x >= 0.05 & Full_results_comb$class == "Proximal",])
G <- D+E
H <- nrow(Full_results_comb[Full_results_comb$P_qvalue.x < 0.05 & Full_results_comb$class == "Distal",])
I <- nrow(Full_results_comb[Full_results_comb$P_qvalue.x >= 0.05 & Full_results_comb$class == "Distal",])
J <- H+I

threebythree <- matrix(ncol = 3, nrow = 6) %>% as.data.frame()
threebythree[,1] <- c(rep("Coding", 2), rep("Proximal", 2), rep("Distal", 2))
threebythree[,2] <- c(A/C,B/C,D/G,E/G,H/J,I/J)
threebythree[,3] <- c("Sig","Non_sig","Sig","Non_sig","Sig","Non_sig")
colnames(threebythree) <- c("Comparison", "Proportion", "Class")

Z <- ggplot(data = threebythree, aes(x=Comparison, y=Proportion, fill=Class)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  ggtitle("Sig vs non-sig CA - within species ")
  ggsave(Z, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Sig_vs_nonSig_CA_change_by_group_within.pdf", width = 10, height = 15)


## Compare variation across time scales
## Estimated accessibility change
df <- Full_results_comb[Full_results_comb$Regulatory_class.x != "Conserved/Ambiguous" & Full_results_comb$Regulatory_class.y != "Conserved/Ambiguous",]
df <- df[,c("P_est.mean.x", "P_est.mean.y")] %>% melt()
df$value <- abs(df$value)
N <- ggplot(df, aes(y=value, x=variable, fill=variable)) +
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
stat_compare_means(method="wilcox.test") +
theme_main() +
xlab("") +
ylab("Estimated CA change") +
scale_x_discrete(labels=c("Within species", "Between species")) +
scale_fill_discrete(guide=FALSE) +
ggtitle("Estimated CA change within vs between species")
ggsave(N, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Est_CA_change_ALL_within_vs_between.pdf", width = 10, height = 15)



## Ratio of sig to non sig
Full_results_comb <- Full_results_comb[,c("P_qvalue.x", "P_qvalue.y")] %>% melt()
A <- nrow(Full_results_comb[Full_results_comb$value < 0.05 & Full_results_comb$variable == "P_qvalue.x",])
B <- nrow(Full_results_comb[Full_results_comb$value >= 0.05 & Full_results_comb$variable == "P_qvalue.x",])
C <- A+B
D <- nrow(Full_results_comb[Full_results_comb$value < 0.05 & Full_results_comb$variable == "P_qvalue.y",])
E <- nrow(Full_results_comb[Full_results_comb$value >= 0.05 & Full_results_comb$variable == "P_qvalue.y",])
G <- D + E

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("Mel_mel", 2), rep("Mel_sim", 2))
twobytwo[,2] <- c(A/C,B/C,D/G,E/G)
twobytwo[,3] <- c("Sig","Non_sig","Sig","Non_sig")
colnames(twobytwo) <- c("Comparison", "Proportion", "Cross")

Y <- ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Cross)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  scale_x_discrete(labels=c("Within species", "Between species")) +
  ggtitle("Sig:non-sig CA - within vs between species")

ggsave(Y, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Sig_vs_nonSig_CA_change_within_vs_between.pdf", width = 11, height = 15)

# Chi Square w/ Yates correction = 0.005

# Same as above but with each class
df <- Full_results_comb_coding[,c("P_est.mean.x", "P_est.mean.y")] %>% melt()
df$value <- abs(df$value)
O <- ggplot(df, aes(y=value, x=variable, fill=variable)) +
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
stat_compare_means(method="wilcox.test") +
theme_main() +
xlab("") +
ylab("Estimated CA change") +
scale_x_discrete(labels=c("Within species", "Between species")) +
scale_fill_discrete(guide=FALSE) +
ggtitle("Estimated CA change within vs between species - CODING")
ggsave(O, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Est_CA_change_CODING_within_vs_between.pdf", width = 10, height = 15)


Full_results_comb_coding <- Full_results_comb_coding[,c("P_qvalue.x", "P_qvalue.y")] %>% melt()
A <- nrow(Full_results_comb_coding[Full_results_comb_coding$value < 0.05 & Full_results_comb_coding$variable == "P_qvalue.x",])
B <- nrow(Full_results_comb_coding[Full_results_comb_coding$value >= 0.05 & Full_results_comb_coding$variable == "P_qvalue.x",])
C <- A+B
D <- nrow(Full_results_comb_coding[Full_results_comb_coding$value < 0.05 & Full_results_comb_coding$variable == "P_qvalue.y",])
E <- nrow(Full_results_comb_coding[Full_results_comb_coding$value >= 0.05 & Full_results_comb_coding$variable == "P_qvalue.y",])
G <- D + E

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("Mel_mel", 2), rep("Mel_sim", 2))
twobytwo[,2] <- c(A/C,B/C,D/G,E/G)
twobytwo[,3] <- c("Sig","Non_sig","Sig","Non_sig")
colnames(twobytwo) <- c("Comparison", "Proportion", "Cross")

Y <- ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Cross)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  scale_x_discrete(labels=c("Within species", "Between species")) +
  ggtitle("Sig:non-sig CA - within vs between species - CODING")

ggsave(Y, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Sig_vs_nonSig_CA_change_within_vs_between_CODING.pdf", width = 11, height = 15)


df <- Full_results_comb_non_coding_dist[,c("P_est.mean.x", "P_est.mean.y")] %>% melt()
df$value <- abs(df$value)
O <- ggplot(df, aes(y=value, x=variable, fill=variable)) +
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
stat_compare_means(method="wilcox.test") +
theme_main() +
xlab("") +
ylab("Estimated CA change") +
scale_x_discrete(labels=c("Within species", "Between species")) +
scale_fill_discrete(guide=FALSE) +
ggtitle("Estimated CA change within vs between species - CODING")
ggsave(O, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Est_CA_change_nonCODING_distal_within_vs_between.pdf", width = 10, height = 15)


Full_results_comb_non_coding_dist <- Full_results_comb_non_coding_dist[,c("P_qvalue.x", "P_qvalue.y")] %>% melt()
A <- nrow(Full_results_comb_non_coding_dist[Full_results_comb_non_coding_dist$value < 0.05 & Full_results_comb_non_coding_dist$variable == "P_qvalue.x",])
B <- nrow(Full_results_comb_non_coding_dist[Full_results_comb_non_coding_dist$value >= 0.05 & Full_results_comb_non_coding_dist$variable == "P_qvalue.x",])
C <- A+B
D <- nrow(Full_results_comb_non_coding_dist[Full_results_comb_non_coding_dist$value < 0.05 & Full_results_comb_non_coding_dist$variable == "P_qvalue.y",])
E <- nrow(Full_results_comb_non_coding_dist[Full_results_comb_non_coding_dist$value >= 0.05 & Full_results_comb_non_coding_dist$variable == "P_qvalue.y",])
G <- D + E

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("Mel_mel", 2), rep("Mel_sim", 2))
twobytwo[,2] <- c(A/C,B/C,D/G,E/G)
twobytwo[,3] <- c("Sig","Non_sig","Sig","Non_sig")
colnames(twobytwo) <- c("Comparison", "Proportion", "Cross")

Y <- ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Cross)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  scale_x_discrete(labels=c("Within species", "Between species")) +
  ggtitle("Sig:non-sig CA - within vs between species - noncoding DISTAL")

ggsave(Y, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Sig_vs_nonSig_CA_change_within_vs_between_nonCODING_dist.pdf", width = 11, height = 15)



df <- Full_results_comb_non_coding_prox[,c("P_est.mean.x", "P_est.mean.y")] %>% melt()
df$value <- abs(df$value)
P <- ggplot(df, aes(y=value, x=variable, fill=variable)) +
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
stat_compare_means(method="wilcox.test") +
theme_main() +
xlab("") +
ylab("Estimated CA change") +
scale_x_discrete(labels=c("Within species", "Between species")) +
scale_fill_discrete(guide=FALSE) +
ggtitle("Estimated CA change within vs between species - noncoding_prox") +
ylim(0,1.5)
ggsave(P, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Est_CA_change_nonCODING_prox_within_vs_between.pdf", width = 10, height = 15)


Full_results_comb_non_coding_prox <- Full_results_comb_non_coding_prox[,c("P_qvalue.x", "P_qvalue.y")] %>% melt()
A <- nrow(Full_results_comb_non_coding_prox[Full_results_comb_non_coding_prox$value < 0.05 & Full_results_comb_non_coding_prox$variable == "P_qvalue.x",])
B <- nrow(Full_results_comb_non_coding_prox[Full_results_comb_non_coding_prox$value >= 0.05 & Full_results_comb_non_coding_prox$variable == "P_qvalue.x",])
C <- A+B
D <- nrow(Full_results_comb_non_coding_prox[Full_results_comb_non_coding_prox$value < 0.05 & Full_results_comb_non_coding_prox$variable == "P_qvalue.y",])
E <- nrow(Full_results_comb_non_coding_prox[Full_results_comb_non_coding_prox$value >= 0.05 & Full_results_comb_non_coding_prox$variable == "P_qvalue.y",])
G <- D + E

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("Mel_mel", 2), rep("Mel_sim", 2))
twobytwo[,2] <- c(A/C,B/C,D/G,E/G)
twobytwo[,3] <- c("Sig","Non_sig","Sig","Non_sig")
colnames(twobytwo) <- c("Comparison", "Proportion", "Cross")

Y <- ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Cross)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  scale_x_discrete(labels=c("Within species", "Between species")) +
  ggtitle("Sig:non-sig CA - within vs between species - noncoding PROX")

ggsave(Y, file="/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/Sig_vs_nonSig_CA_change_within_vs_between_nonCODING_prox.pdf", width = 11, height = 15)


## WHAT IS THE MODE OF ACCESSIBILITY DIVERGENCE ACROSS CLASSES AND TIME?
# Plot % cis across groups
A <- Full_results_comb[Full_results_comb$Regulatory_class.x != "Conserved/Ambiguous",] %>%
ggplot(aes(y=perc_cis.x, x=class, fill=class)) +
geom_boxplot(notch=TRUE) +
xlab("") +
ylab("% cis") +
scale_fill_discrete(guide=FALSE) +
theme_main() +
ggtitle("Within species comparison of % cis across groups")
ggsave(A, file = "~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/perc_cis_within_ATAC_across_groups.pdf", width = 10, height = 15)

B <- Full_results_comb[Full_results_comb$Regulatory_class.y != "Conserved/Ambiguous",] %>%
ggplot(aes(y=perc_cis.y, x=class, fill=class)) +
geom_boxplot(notch=TRUE) +
xlab("") +
ylab("% cis") +
scale_fill_discrete(guide=FALSE) +
theme_main() +
ggtitle("Between species comparison of % cis across groups")
ggsave(B, file = "~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/perc_cis_between_ATAC_across_groups.pdf", width = 10, height = 15)


stat.test <- Full_results_comb[Full_results_comb$Regulatory_class.x != "Conserved/Ambiguous",] %>%
  wilcox_test(perc_cis.x ~ class)

stat.test <- Full_results_comb[Full_results_comb$Regulatory_class.y != "Conserved/Ambiguous",] %>%
  wilcox_test(perc_cis.y ~ class)

# Plot % cis across time for coding vs non-coding
## non-coding
Full_results_comb_div <- Full_results_comb_non_coding[Full_results_comb_non_coding$Regulatory_class.x != "Conserved/Ambiguous" & Full_results_comb_non_coding$Regulatory_class.y != "Conserved/Ambiguous",]
Perc_cis <- cbind(Full_results_comb_div $perc_cis.x, Full_results_comb_div $perc_cis.y)
colnames(Perc_cis) <- c("Mel_Mel", "Mel_Sim")
Perc_cis_melt <- melt(Perc_cis)

R  <- ggplot(Perc_cis_melt, aes(x=Var2, y=value, fill=Var2)) +
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
scale_x_discrete(labels=c("Within species","Between species")) +
xlab("") +
ylab("% cis") +
scale_fill_discrete(guide=FALSE) +
theme_main() +
stat_compare_means(method = "wilcox.test")
ggsave(R, file = "~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/perc_cis_within_and_between_ATAC_non_coding.pdf", width = 10, height = 15)

## coding
Full_results_comb_div <- Full_results_comb_coding[Full_results_comb_coding$Regulatory_class.x != "Conserved/Ambiguous" & Full_results_comb_coding$Regulatory_class.y != "Conserved/Ambiguous",]
Perc_cis <- cbind(Full_results_comb_div $perc_cis.x, Full_results_comb_div $perc_cis.y)
colnames(Perc_cis) <- c("Mel_Mel", "Mel_Sim")
Perc_cis_melt <- melt(Perc_cis)

R  <- ggplot(Perc_cis_melt, aes(x=Var2, y=value, fill=Var2)) +
geom_violin() +
geom_boxplot(notch=TRUE, width = 0.2) +
scale_x_discrete(labels=c("Within species","Between species")) +
xlab("") +
ylab("% cis") +
scale_fill_discrete(guide=FALSE) +
theme_main() +
stat_compare_means(method = "wilcox.test")
ggsave(R, file = "~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/perc_cis_within_and_between_ATAC_coding.pdf", width = 10, height = 15)

## HOW DOES CIS&TRANS DIFFER ACROSS TIME AND SPACE?

## All regions
df <- Full_results_comb_div
x <- nrow(df[df$Direction.x == "Opposing",])
y <- nrow(df[df$Direction.x == "Reinforcing",])
intra_opposing <- x/(x+y)
intra_reinforcing <- y/(x+y)

a <- nrow(df[df$Direction.y == "Opposing",])
b <- nrow(df[df$Direction.y == "Reinforcing",])
inter_opposing <- a/(a+b)
inter_reinforcing <- b/(a+b)

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("Mel-Mel", 2), rep("Mel-Sim", 2))
twobytwo[,2] <- c(intra_opposing,intra_reinforcing,inter_opposing,inter_reinforcing)
twobytwo[,3] <- c("Cis-trans opposing","Cis-trans reinforcing","Cis-trans opposing","Cis-trans reinforcing")
colnames(twobytwo) <- c("Comparison", "Proportion", "Cis_trans_Directionality")

B <- ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Cis_trans_Directionality)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  scale_x_discrete(labels=c("Within","Between")) +
  ggtitle("Cis-trans oppo vs reinforce within & between species")

ggsave(B, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/cis_trans_opposing_reinforcing_ATAC.pdf", width = 10, height = 15)

## Coding regions
Full_results_comb_div <- Full_results_comb_coding[Full_results_comb_coding$Regulatory_class.x != "Conserved/Ambiguous" & Full_results_comb_coding$Regulatory_class.y != "Conserved/Ambiguous",]
df <- Full_results_comb_div
x <- nrow(df[df$Direction.x == "Opposing",])
y <- nrow(df[df$Direction.x == "Reinforcing",])
intra_opposing <- x/(x+y)
intra_reinforcing <- y/(x+y)

a <- nrow(df[df$Direction.y == "Opposing",])
b <- nrow(df[df$Direction.y == "Reinforcing",])
inter_opposing <- a/(a+b)
inter_reinforcing <- b/(a+b)

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("Mel-Mel", 2), rep("Mel-Sim", 2))
twobytwo[,2] <- c(intra_opposing,intra_reinforcing,inter_opposing,inter_reinforcing)
twobytwo[,3] <- c("Cis-trans opposing","Cis-trans reinforcing","Cis-trans opposing","Cis-trans reinforcing")
colnames(twobytwo) <- c("Comparison", "Proportion", "Cis_trans_Directionality")

B <- ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Cis_trans_Directionality)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  scale_x_discrete(labels=c("Within","Between")) +
  ggtitle("Cis-trans oppo vs reinforce within & between species")

ggsave(B, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/cis_trans_opposing_reinforcing_ATAC_coding.pdf", width = 10, height = 15)

## Distal
Full_results_comb_div <- Full_results_comb_non_coding_dist[Full_results_comb_non_coding_dist$Regulatory_class.x != "Conserved/Ambiguous" & Full_results_comb_non_coding_dist$Regulatory_class.y != "Conserved/Ambiguous",]
df <- Full_results_comb_div
x <- nrow(df[df$Direction.x == "Opposing",])
y <- nrow(df[df$Direction.x == "Reinforcing",])
intra_opposing <- x/(x+y)
intra_reinforcing <- y/(x+y)

a <- nrow(df[df$Direction.y == "Opposing",])
b <- nrow(df[df$Direction.y == "Reinforcing",])
inter_opposing <- a/(a+b)
inter_reinforcing <- b/(a+b)

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("Mel-Mel", 2), rep("Mel-Sim", 2))
twobytwo[,2] <- c(intra_opposing,intra_reinforcing,inter_opposing,inter_reinforcing)
twobytwo[,3] <- c("Cis-trans opposing","Cis-trans reinforcing","Cis-trans opposing","Cis-trans reinforcing")
colnames(twobytwo) <- c("Comparison", "Proportion", "Cis_trans_Directionality")

B <- ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Cis_trans_Directionality)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  scale_x_discrete(labels=c("Within","Between")) +
  ggtitle("Cis-trans oppo vs reinforce within & between species")

ggsave(B, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/cis_trans_opposing_reinforcing_ATAC_non_coding_dist.pdf", width = 10, height = 15)

## Proximal
Full_results_comb_div <- Full_results_comb_non_coding_prox[Full_results_comb_non_coding_prox$Regulatory_class.x != "Conserved/Ambiguous" & Full_results_comb_non_coding_prox$Regulatory_class.y != "Conserved/Ambiguous",]
df <- Full_results_comb_div
x <- nrow(df[df$Direction.x == "Opposing",])
y <- nrow(df[df$Direction.x == "Reinforcing",])
intra_opposing <- x/(x+y)
intra_reinforcing <- y/(x+y)

a <- nrow(df[df$Direction.y == "Opposing",])
b <- nrow(df[df$Direction.y == "Reinforcing",])
inter_opposing <- a/(a+b)
inter_reinforcing <- b/(a+b)

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("Mel-Mel", 2), rep("Mel-Sim", 2))
twobytwo[,2] <- c(intra_opposing,intra_reinforcing,inter_opposing,inter_reinforcing)
twobytwo[,3] <- c("Cis-trans opposing","Cis-trans reinforcing","Cis-trans opposing","Cis-trans reinforcing")
colnames(twobytwo) <- c("Comparison", "Proportion", "Cis_trans_Directionality")

B <- ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Cis_trans_Directionality)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  scale_x_discrete(labels=c("Within","Between")) +
  ggtitle("Cis-trans oppo vs reinforce within & between species")

ggsave(B, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/cis_trans_opposing_reinforcing_ATAC_non_coding_prox.pdf", width = 10, height = 15)
