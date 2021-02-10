# originally writing this to analyze ATAC-seq data, but could be re-used for other data types
### question: is there a correlation between chromatin accessibility changes by proximity?
## approach: read in differentially accessible peaks and look at scatter plot between focal peak vs "peak to the right"

#########################
##Set master plot theme##
#########################
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

# load libraries and functions
library(stringr)
library(dplyr)
library(plyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(rstatix)

# read in datasets
diff_peaks <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATACseq/Diffbind_file/ZHR_Z30_diff_accessible_final_no_quotes_sort.bed", header = F)
diff_peaks_clean <- diff_peaks[,c(1,2,3,9)]
colnames(diff_peaks_clean) <- c("chrom", "start", "end", "fold")

# make dummy closest set (had to play around in and brute force it with nano) and append
peaks_down_one <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATACseq/Diffbind_file/ZHR_Z30_diff_accessible_final_no_quotes_sort_shifted_down1.bed", header = F)
temp_col <- peaks_down_one[,c(1,2,3,9)]
colnames(temp_col) <- c("chrom_closest", "start_closest", "end_closest", "fold_closest")
final_data_frame <- as.data.frame(cbind(diff_peaks_clean, temp_col))

# check distance between peaks
final_data_frame$diff <- (final_data_frame$start_closest - final_data_frame$start)

# eliminante weird points that are negative difference (e.g. chromosome ends) & outliers by force *look at later
final_data_frame_subset <- as.data.frame(subset(final_data_frame_subset, diff > 0 & fold > 0 & fold < 2 & fold_closest > 0 & fold_closest < 2))

# plot to see all data
ggplot(final_data_frame_subset, aes(x=abs(fold), y=abs(fold_closest))) +
  geom_point()

ggplot(final_data_frame_subset, aes(x=diff)) +
  geom_density() +
  xlim(0, 30000)

# analyze correlations by distance from peaks
## overlapping sets:
l <- seq(from = 2000, to = 30000, by = 500)
index <- c(1:length(l))
l_cor <- rep(0,length(l))

for (k in index) {
  df <- subset(final_data_frame_subset, final_data_frame_subset$diff > 0 & final_data_frame_subset$diff < l[k])
  l_cor[k] <- cor(df$fold,df$fold_closest,method="pearson")
}

### make data tables and plot
cor_by_dist <- NULL
cor_by_dist$cor <- l_cor
cor_by_dist$distance <- l/1000
cor_by_dist <- as.data.frame(cor_by_dist)

cor_by_dist_plot <- ggplot(cor_by_dist, aes(x=distance, y=cor)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ylab("Pearson cor: focal vs nearest peak fold change") +
  scale_x_continuous(name="Distance from focal peak (Kb)", breaks=(seq(2,30,2))) +
  theme_main()

ggsave(cor_by_dist_plot, file = "~/Documents/Wittkopp_lab/AS_ATACseq/Figures/ZHR_Z30_pilot_cor_by_dist_overlapping.pdf", width = 15, height = 15)


## non-overlapping sets:
l_min <- seq(from = 0, to = 9000, by = 1000)
l_max <- seq(from = 1000, to = 10000, by = 1000)
index_2 <- c(1:length(l_min))
l_cor_2 <- rep(0,length(l_min))

for (k in index_2) {
  df_2 <- subset(final_data_frame_subset, final_data_frame_subset$diff > l_min[k] & final_data_frame_subset$diff < l_max[k])
  l_cor_2[k] <- cor(df_2$fold,df_2$fold_closest,method="pearson")
}

cor_by_dist_2 <- NULL
cor_by_dist_2$cor <- l_cor_2
cor_by_dist_2$distance <- l_max/1000
cor_by_dist_2 <- as.data.frame(cor_by_dist_2)

cor_by_dist_plot_2 <- ggplot(cor_by_dist_2, aes(x=distance, y=cor)) +
  geom_point() +
  stat_smooth(method = "lm") +
  ylab("Pearson cor: focal vs nearest peak fold change") +
  scale_x_continuous(name="Distance from focal peak (Kb)", breaks=(seq(1,10,1))) +
  theme_main()

ggsave(cor_by_dist_plot_2, file = "~/Documents/Wittkopp_lab/AS_ATACseq/Figures/ZHR_Z30_pilot_cor_by_dist_nonoverlapping.pdf", width = 15, height = 15)

cor_2 <- NULL
cor_2 <- as.data.frame(subset(final_data_frame_subset, final_data_frame_subset$diff > 19500 & final_data_frame_subset$diff < 20000))

ggplot(cor_2, aes(x=fold, y=fold_closest)) +
  geom_point()






######
library(plyr)
library(dplyr)
library(ggplot2)

df <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/ZHR_Z30_TSIM_ATAC_20min_peaks_SNPs_intersect.bed", header = F)
df <- df[,c(1:4)]
colnames(df) <- c("chrom", "start", "end", "paste_locus")
df <- ddply(df,.(chrom,start,end,paste_locus),nrow)
colnames(df) <- c("chrom", "start", "end", "Paste_locus", "SNP_num")
df$length <- df$end - df$start
df$SNP_length_norm <- df$SNP_num/df$length * 1000
df_final <- df[,c(4,7)]

# merge with data frame

data <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/ZHR_Z30_TSIM_universal_peakset_outputs/Full_results_output_ZHR_Z30_TSIM_ATAC_20min.txt", header = T)
data_new <- full_join(df_final, data, by = "Paste_locus")

A <-
data_new[data_new$Regulatory_class.x == "Cis" | data_new$Regulatory_class.x == "Trans" | data_new$Regulatory_class.x == "Conserved/Ambiguous",] %>%
ggplot(aes(x=Regulatory_class.x, fill=Regulatory_class.x, y=SNP_length_norm)) +
geom_boxplot(with=0.5, notch = TRUE) +
theme_main() +
ylab("Region variant density (variant # / 1 kB)") +
xlab("") +
scale_fill_discrete(guide=FALSE) +
scale_x_discrete(labels=c("Cis","Conserved", "Trans"))

ggsave(A, file = "~/Documents/Wittkopp_lab/AS_ATACseq/Figures/ZHR_Z30_variant_reg_class.pdf", width = 10, height = 15)


stat.test <- data_new[data_new$Regulatory_class.x == "Cis" | data_new$Regulatory_class.x == "Trans" | data_new$Regulatory_class.x == "Conserved/Ambiguous",] %>%
  wilcox_test(Regulatory_class.x ~ SNP_length_norm) %>%
  add_significance()
stat.test


Y <-
data_new[data_new$Regulatory_class.x != "Conserved/Ambiguous",] %>%
ggplot() +
geom_point(aes(x=SNP_length_norm, y=perc_cis.x)) +
theme_main() +
xlab("Region variant density (variant # / 1 kB)") +
ylab("% Cis")
ggsave(Y, file = "~/Documents/Wittkopp_lab/AS_ATACseq/Figures/ZHR_Z30_variant_perc_cis.pdf", width = 15, height = 15)





data_new[data_new$Regulatory_class.x == "Cis",] %>%
ggplot() +
geom_point(aes(x=SNP_length_norm, y=abs(P_est.mean.x)))


data_new

## Assign inDEGREE quartile to each gene

data_new$perc_cis.x_rank <- "NA"

for (i in 1:nrow(data_new)) {

if (data_new$perc_cis.x[i] <= 25){

	data_new$perc_cis.x_rank[i] <- "0_25"

} else if (data_new$perc_cis.x[i] > 25 & data_new$perc_cis.x[i] <= 50){

	data_new$perc_cis.x_rank[i] <- "25_50"

} else if (data_new$perc_cis.x[i] > 50 & data_new$perc_cis.x[i] <= 75){

	data_new$perc_cis.x_rank[i] <- "50_75"

} else if (data_new$perc_cis.x[i] > 75){

	data_new$perc_cis.x_rank[i] <- "75_100"

}
}



A <- ggplot(data_new, aes(x=perc_cis.x_rank, fill=perc_cis.x_rank, y=SNP_length_norm)) +
geom_boxplot()
