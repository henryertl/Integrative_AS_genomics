library(ggplot2)
library(plyr)
library(dplyr)

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


# read in file
ZHR_Z30_regions_ALL_SNPs_GRh_final <- read.delim("./AS_ATAC_RNA_2020_10_1/Grh_data/ZHR_Z30_GRH_SNPs_count_FINALALL.txt", header = T)

# plot cis-trans to set boundaries
cis_trans_ATAC_CPM <- ZHR_Z30_regions_ALL_SNPs_GRh_final  %>%

ggplot(aes(x = P_est.mean, y = H_est.mean)) +
geom_point(alpha = 0.3) +
geom_bin2d(bins = 500) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + theme_main() + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-2, 2) + ylim(-2, 2)

cis_trans_ATAC_CPM_chrom <- ZHR_Z30_regions_ALL_SNPs_GRh_final  %>%

ggplot(aes(x = P_est.mean, y = H_est.mean)) +
geom_point(alpha = 0.3) +
geom_bin2d(bins = 100) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + theme_main() + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-2, 2) + ylim(-2, 2) +
facet_wrap(~chrom)

## chrom2L is our issue
ZHR_Z30_regions_ALL_SNPs_GRh_final %>%
ggplot(aes(x=H_est.mean, color = chrom, fill = chrom)) +
geom_density()

# look at divergence per chrom
E <- ZHR_Z30_regions_ALL_SNPs_GRh_final %>%

ggplot(aes(x=chrom, y=abs(P_est.mean), fill=chrom))+
geom_boxplot(notch=T) +
ylim(0,1) +
theme_main() +
ylab("Estimted accessibility divergence") +
xlab("") +
scale_fill_discrete(guide=FALSE)


F <- ZHR_Z30_regions_ALL_SNPs_GRh_final %>%

ggplot(aes(x=chrom, y=abs(total_region_snps), fill=chrom))+
geom_boxplot(notch=T) +
theme_main() +
ylab("# of nucleotide variants in region") +
xlab("") +
scale_fill_discrete(guide=FALSE) +


ZHR_Z30_regions_ALL_SNPs_GRh_final[ZHR_Z30_regions_ALL_SNPs_GRh_final$chrom != "chr4",] %>%
ggplot(aes(x=end, fill = chrom)) +
geom_density() +
facet_wrap(~chrom) +
scale_fill_discrete(guide=FALSE)


# is chr2L more seq divergent overall (outside of CA regions?)

snps <- read.delim("/Users/wittkopp_member/Code/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/zhr_z30_SNPs_dm6_genome.bed", header = T)
colnames(snps) <- c("chrom", "start", "end")
snps_chrom <- snps$chrom %>% as.data.frame()
colnames(snps_chrom) <- "chrom"

snps_chrom_count <- ddply(snps_chrom,.(chrom),nrow)
snps_chrom_count$chrom_size <- c(23513712, 25286936, 28110227, 32079331, 1348131, 23542271)
colnames(snps_chrom_count) <- c("chrom", "snps", "chrom_size")

snps_chrom_count$ratio <- snps_chrom_count$snps/snps_chrom_count$chrom_size

snps_chrom_count %>%

ggplot(aes(x=chrom, y=ratio)) +
geom_col() +
ylab("snps/chr size")
