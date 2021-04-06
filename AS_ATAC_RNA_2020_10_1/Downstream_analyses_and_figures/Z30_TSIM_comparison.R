setwd("/Users/henryertl/Documents/Devs/Integrative_AS_genomics")
setwd("/Users/wittkopp_member/Code")


# read in both Full_results_output fiels

Z30 <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.txt", header = T) %>% na.omit() %>% unique()
TSIM <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_TSIM_ATAC_20min_downsamp_overlap.txt", header = T) %>% na.omit() %>% unique()


colnames(Z30) <- paste(colnames(Z30), "Z30", sep = "_")
colnames(Z30)[4] <- "Paste_locus"
colnames(TSIM) <- paste(colnames(TSIM), "TSIM", sep = "_")
colnames(TSIM)[4] <- "Paste_locus"

Z30_TSIM <- join_all(list(Z30, TSIM), type = "full", by = "Paste_locus") %>% na.omit() %>% unique()

Z30_TSIM$class_Z30 <- factor(Z30_TSIM$class_Z30,levels = c("start", "end", "inter", "intra"))

write.table(Z30_TSIM, file = "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/ZHR_Z30_TSIM_Full_results_output_ATAC.txt", quote = F, sep = "\t", row.names = F)


# plot percent cis for each class within vs between species
A <- Z30_TSIM[(Z30_TSIM$P_qvalue_TSIM < 0.05 | Z30_TSIM$H_qvalue_TSIM < 0.05) & (Z30_TSIM$P_qvalue_Z30 < 0.05 | Z30_TSIM$H_qvalue_Z30 < 0.05),] %>%
melt(id.vars = "class_Z30", measure.vars = c("perc_cis_Z30", "perc_cis_TSIM")) %>%
ggplot(aes(x=variable, y=value, fill=variable)) +
geom_boxplot(notch=TRUE) +
theme_main() +
ylab("CA divergence due to cis changes (percent)") +
xlab("") +
scale_x_discrete(labels=c("Within\nspecies", "Between\nspecies")) +
scale_fill_discrete(guide=F) +
facet_wrap(~class_Z30, nrow=2, , labeller = labeller(class_Z30=
    c("start" = "txStart",
      "end" = "txEnd",
      "inter" = "intergenic",
      "intra" = "intragenic")
  ))

ggsave(A, file = "./AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/perc_cis_withinVSbetween.pdf", width = 4.5)

C <- Z30_TSIM[(Z30_TSIM$P_qvalue_TSIM < 0.05 | Z30_TSIM$H_qvalue_TSIM < 0.05) & (Z30_TSIM$P_qvalue_Z30 < 0.05 | Z30_TSIM$H_qvalue_Z30 < 0.05),] %>%
melt(measure.vars = c("perc_cis_Z30", "perc_cis_TSIM")) %>%
ggplot(aes(x=variable, y=value, fill=variable)) +
geom_boxplot(notch=TRUE) +
theme_main() +
ylab("CA divergence due to cis changes (percent)") +
xlab("") +
scale_x_discrete(labels=c("Within\nspecies", "Between\nspecies")) +
scale_fill_discrete(guide=F)

ggsave(C, file = "./AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/perc_cis_withinVSbetween_ALL.pdf", width = 3)


B <- Z30_TSIM %>%
melt(id.vars = "class_Z30", measure.vars = c("P_est.mean_Z30", "P_est.mean_TSIM")) %>%
ggplot(aes(x=variable, y=abs(value), fill=variable)) +
geom_boxplot(notch=TRUE) +
ylim(0,1) +
theme_main() +
ylab("Estimated CA divergence ") +
xlab("") +
scale_x_discrete(labels=c("Within\nspecies", "Between\nspecies")) +
scale_fill_discrete(guide=F) +
facet_wrap(~class_Z30, nrow=2, , labeller = labeller(class_Z30=
    c("start" = "txStart",
      "end" = "txEnd",
      "inter" = "intergenic",
      "intra" = "intragenic")
  ))

ggsave(B, file = "./AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/Est_CA_divergence_withinVSbetween.pdf", width = 4.5)

D <- Z30_TSIM %>%
melt(measure.vars = c("H_est.mean_Z30", "H_est.mean_TSIM")) %>%
ggplot(aes(x=variable, y=abs(value), fill=variable)) +
geom_boxplot(notch=TRUE) +
ylim(0,1) +
theme_main() +
ylab("Estimated CA divergence ") +
xlab("") +
scale_x_discrete(labels=c("Within\nspecies", "Between\nspecies")) +
scale_fill_discrete(guide=F)

ggsave(D, file = "./AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/Est_CA_divergence_withinVSbetween_ALL.pdf", width = 3)
