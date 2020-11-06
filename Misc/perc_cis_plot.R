## from Bayes outputs of ZHR - TSIM and ZHR - Z30 comparisons, compare percent cis

# RNA_seq

TSIM <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_TSIM_RNA_20min.txt", header = T)
Z30 <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_1000max.txt", header = T)

Full_results_comb <- merge(Z30, TSIM, by = "Paste_locus")
Full_results_comb <- Full_results_comb[-c(1)]
write.table(Full_results_comb, file="~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/ZHR_Z30_TSIM_Full_results_output.txt", sep = "\t", row.names = F, quote = F)

Perc_cis <- cbind(Full_results_comb$perc_cis.x, Full_results_comb$perc_cis.y)
colnames(Perc_cis) <- c("Mel_Mel", "Mel_Sim")
Perc_cis_melt <- melt(Perc_cis)

stat.test <- Perc_cis_melt  %>%
  wilcox_test(value ~ Var2, paired = TRUE) %>%
  add_significance()
stat.test

RNA_perc_cis <- ggplot(Perc_cis_melt, aes(x=Var2, y=value, fill=Var2)) +
  geom_violin() +
  xlab("") +
  ylab("% cis") +
  #stat_summary(fun.y='mean', geom='point', size=2, col='red') +
  geom_boxplot(width=0.1) +
  theme_main() +
  theme(legend.position="none") +
  ggtitle("Gene Expression")

# ATAC_seq

TSIM <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/ZHR_Z30_TSIM_universal_peakset_outputs/Full_results_output_ZHR_TSIM_ATAC_20min.txt", header = T)
Z30 <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/ZHR_Z30_TSIM_universal_peakset_outputs/Full_results_output_ZHR_Z30_ATAC_20min.txt", header = T)

Full_results_comb <- merge(Z30, TSIM, by = "Paste_locus")
Full_results_comb <- Full_results_comb[-c(1)]
write.table(Full_results_comb, file="~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/ZHR_Z30_TSIM_universal_peakset_outputs/ZHR_Z30_TSIM_Full_results_output.txt", sep = "\t", row.names = F, quote = F)

Perc_cis <- cbind(Full_results_comb$perc_cis.x, Full_results_comb$perc_cis.y)
colnames(Perc_cis) <- c("Mel_Mel", "Mel_Sim")
Perc_cis_melt <- melt(Perc_cis)

stat.test <- Perc_cis_melt  %>%
  wilcox_test(value ~ Var2, paired = TRUE) %>%
  add_significance()
stat.test

ATAC_perc_cis <- ggplot(Perc_cis_melt, aes(x=Var2, y=value, fill=Var2)) +
  geom_violin() +
  xlab("") +
  ylab("% cis") +
  #stat_summary(fun.y='mean', geom='point', size=2, col='red') +
  geom_boxplot(width=0.1) +
  theme_main() +
  theme(legend.position="none") +
  ggtitle("Chromatin Accessibility")

final <- plot_grid(RNA_perc_cis, ATAC_perc_cis)
ggsave(final, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_RNA_comp/perc_cis_ZHR_Z30_TSIM_RNA_ATAC.pdf", width = 15, height = 15)
