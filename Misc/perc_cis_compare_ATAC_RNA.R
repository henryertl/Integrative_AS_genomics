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
  axis.text = element_text(size = 20),
  axis.title = element_text(size = 30),
  strip.text = element_text(size = 20),
  legend.text= element_text(size = 20),
  legend.title = element_text(size = 15),
  plot.title = element_text(size = 17, face = "bold")

)
}

#######################################################################
###### Generate base plots - % cis accumulation and differences #######
#######################################################################

# Read in Bayes output files
ZHR_Z30_TSIM_RNA <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/ZHR_Z30_TSIM_Full_results_output.txt", header = T)
ZHR_Z30_TSIM_RNA$data_type <- "RNA"
ZHR_Z30_TSIM_ATAC <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/ZHR_Z30_TSIM_Full_results_output_ATAC.txt", header = T)
ZHR_Z30_TSIM_ATAC$data_type <- "ATAC"

ZHR_Z30_TSIM_RNA_min <- ZHR_Z30_TSIM_RNA[(ZHR_Z30_TSIM_RNA$P_qvalue.y < 0.05 |ZHR_Z30_TSIM_RNA$H_qvalue.y < 0.05) & (ZHR_Z30_TSIM_RNA$P_qvalue.x < 0.05 | ZHR_Z30_TSIM_RNA$H_qvalue.x < 0.05),]
ZHR_Z30_TSIM_RNA_div <- cbind(ZHR_Z30_TSIM_RNA_min$P_est.mean.x, ZHR_Z30_TSIM_RNA_min$P_est.mean.y) %>% as.data.frame() %>% na.omit()
ZHR_Z30_TSIM_RNA_perc_cis <- cbind(ZHR_Z30_TSIM_RNA_min$perc_cis.x, ZHR_Z30_TSIM_RNA_min$perc_cis.y) %>% as.data.frame() %>% na.omit()
colnames(ZHR_Z30_TSIM_RNA_div) <- c("P_est.mean_Z30", "P_est.mean_TSIM")
colnames(ZHR_Z30_TSIM_RNA_perc_cis) <- c("perc_cis_Z30", "perc_cis_TSIM")


ZHR_Z30_TSIM_ATAC_min <- ZHR_Z30_TSIM_ATAC[(ZHR_Z30_TSIM_ATAC$P_qvalue_TSIM < 0.05 |ZHR_Z30_TSIM_ATAC$H_qvalue_TSIM < 0.05) & (ZHR_Z30_TSIM_ATAC$P_qvalue_Z30 < 0.05 | ZHR_Z30_TSIM_ATAC$H_qvalue_Z30 < 0.05),]
ZHR_Z30_TSIM_ATAC_div <- cbind(ZHR_Z30_TSIM_ATAC_min$P_est.mean_Z30, ZHR_Z30_TSIM_ATAC_min$P_est.mean_TSIM) %>% as.data.frame() %>% na.omit()
ZHR_Z30_TSIM_ATAC_perc_cis <- cbind(ZHR_Z30_TSIM_ATAC_min$perc_cis_Z30, ZHR_Z30_TSIM_ATAC_min$perc_cis_TSIM) %>% as.data.frame() %>% na.omit()
colnames(ZHR_Z30_TSIM_ATAC_div) <- c("P_est.mean_Z30", "P_est.mean_TSIM")
colnames(ZHR_Z30_TSIM_ATAC_perc_cis) <- c("perc_cis_Z30", "perc_cis_TSIM")

ZHR_Z30_TSIM_ATAC_div_melt <- melt(ZHR_Z30_TSIM_ATAC_div)
ZHR_Z30_TSIM_ATAC_div_melt$datatype <- "ATAC"
ZHR_Z30_TSIM_RNA_div_melt <- melt(ZHR_Z30_TSIM_RNA_div)
ZHR_Z30_TSIM_RNA_div_melt$datatype <- "RNA"
ZHR_Z30_TSIM_ATAC_RNA_div <- rbind(ZHR_Z30_TSIM_ATAC_div_melt, ZHR_Z30_TSIM_RNA_div_melt)

ZHR_Z30_TSIM_ATAC_perc_cis_melt <- melt(ZHR_Z30_TSIM_ATAC_perc_cis)
ZHR_Z30_TSIM_ATAC_perc_cis_melt$datatype <- "ATAC"
ZHR_Z30_TSIM_RNA_perc_cis_melt <- melt(ZHR_Z30_TSIM_RNA_perc_cis)
ZHR_Z30_TSIM_RNA_perc_cis_melt$datatype <- "RNA"
ZHR_Z30_TSIM_ATAC_RNA_perc_cis <- rbind(ZHR_Z30_TSIM_ATAC_perc_cis_melt, ZHR_Z30_TSIM_RNA_perc_cis_melt)

## Divergence
R  <- ZHR_Z30_TSIM_ATAC_RNA_div %>%
ggplot(aes(x=variable, y=abs(value), fill=variable)) +
geom_boxplot(notch=TRUE) +
ylim(0,1) +
theme_main() +
ylab("Estimated CA divergence ") +
xlab("") +
scale_x_discrete(labels=c("Within\nspecies", "Between\nspecies")) +
scale_fill_discrete(guide=F) +
facet_wrap(~datatype, nrow=1)
ggsave(R, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/Est_CA_divergence_withinVSbetween_ALL_RNA_ATAC.pdf", width = 5)


S <- ZHR_Z30_TSIM_ATAC_RNA_perc_cis %>%
ggplot(aes(x=variable, y=abs(value), fill=variable)) +
geom_boxplot(notch=TRUE) +
theme_main() +
ylab("perc cis") +
xlab("") +
scale_x_discrete(labels=c("Within\nspecies", "Between\nspecies")) +
scale_fill_discrete(guide=F) +
facet_wrap(~datatype, nrow=1)
ggsave(S, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/perc_cis_withinVSbetween_ALL_RNA_ATAC.pdf", width = 5)
