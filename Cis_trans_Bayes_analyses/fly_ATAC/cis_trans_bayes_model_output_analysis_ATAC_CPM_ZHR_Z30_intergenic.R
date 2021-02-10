#############
##Libraries##
#############
library(dplyr)
library(plyr)
library(cowplot)
library(magrittr)
library(ggplot2)


#########################
##Set master plot theme##
#########################
theme_main <- function() {
  theme_bw() +
  theme(
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 20),
  strip.text = element_text(size = 20),
  legend.text= element_text(size = 15),
  legend.title = element_text(size = 20),
  plot.title = element_text(size = 25, face = "bold")
)
}

#####################
##Read primary data##
#####################
full_dataset <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Data_tables/ZHR_Z30_ATAC_intergenic_CPM_final_dm6_20min.txt", header = T)
colnames(full_dataset) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")
full_dataset$Paste_locus <- paste(full_dataset$chrom, full_dataset$start, full_dataset$end, sep = "_")

Parental_data <- full_dataset[, c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3")]

Hybrid_data <- full_dataset[, c("chrom", "start", "end", "Paste_locus", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")]


####################
##Combine datasets##
####################
Parental_results <- read.table("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/Parental_test_output_ZHR_Z30_ATAC_CPM_macs2_20min_intergenic.txt", header = T)
Hybrid_results <- read.table("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/Hybrid_test_output_ZHR_Z30_ATAC_CPM_macs2_20min_intergenic.txt", header = T)
Parental_hybrid_results <- read.table("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/Parental_Hybrid_test_output_ZHR_Z30_ATAC_CPM_macs2_20min_intergenic.txt", header = T)

Full_results_output <- join_all(list(Parental_data, Hybrid_data, Parental_results, Hybrid_results, Parental_hybrid_results), by = 'Paste_locus', type = 'full')

Full_results_output <- na.omit(Full_results_output)

##Apply FDR correction
##First plot distribution of p-values to check for any weirdness
Parent_P_plot <- ggplot(Full_results_output, aes(x = P_p_value)) + geom_histogram(bins = 100) + ggtitle("Parents")
Hybrid_P_plot <- ggplot(Full_results_output, aes(x = H_p_value)) + geom_histogram(bins = 100) + ggtitle("Hybrids")
Parent_hybrid_P_plot <- ggplot(Full_results_output, aes(x = H_P_p_value)) + geom_histogram(bins = 100) + ggtitle("Hybrid-Parents")

p_vals <- plot_grid(Parent_P_plot, Hybrid_P_plot, Parent_hybrid_P_plot, nrow = 3)
ggsave(p_vals, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/p_vals_ZHR_Z30_sub_ATAC_20min)intergenic.pdf", width = 15, height = 15)


##If all is well, run FDR correction
Full_results_output$P_qvalue <- p.adjust(Full_results_output$P_p_value, method = "BY")
Full_results_output$H_qvalue <- p.adjust(Full_results_output$H_p_value, method = "BY")
Full_results_output$P_H_qvalue <- p.adjust(Full_results_output$H_P_p_value, method = "BH")

####################################
##Get consistent allele directions##
####################################
Full_results_output$Direction_parent <- NA

for (i in 1:nrow(Full_results_output)) {

if (Full_results_output$P1_1[i] > Full_results_output$P2_1[i] & Full_results_output$P1_2[i] > Full_results_output$P2_2[i] & Full_results_output$P1_3[i] > Full_results_output$P2_3[i]){

	Full_results_output$Direction_parent[i] <- "P1"

} else if (Full_results_output$P1_1[i] < Full_results_output$P2_1[i] & Full_results_output$P1_2[i] < Full_results_output$P2_2[i] & Full_results_output$P1_3[i] < Full_results_output$P2_3[i]){

	Full_results_output$Direction_parent[i] <- "P2"

} else {Full_results_output$Direction_parent[i] <- "Ambig"}
}


Full_results_output$Direction_hybrid <- NA

for (i in 1:nrow(Full_results_output)) {

if (Full_results_output$HYB_1_P1[i] > Full_results_output$HYB_1_P2[i] & Full_results_output$HYB_2_P1[i] > Full_results_output$HYB_2_P2[i] & Full_results_output$HYB_3_P1[i] > Full_results_output$HYB_3_P2[i]){

	Full_results_output$Direction_hybrid[i] <- "P1"

} else if (Full_results_output$HYB_1_P1[i] < Full_results_output$HYB_1_P2[i] & Full_results_output$HYB_2_P1[i] < Full_results_output$HYB_2_P2[i] & Full_results_output$HYB_3_P1[i] < Full_results_output$HYB_3_P2[i]){

	Full_results_output$Direction_hybrid[i] <- "P2"

} else {Full_results_output$Direction_hybrid[i] <- "Ambig"}
}

##########################################################################
### Get Parental / Hybrid ratio to classify opposing vs same cis+trans ###
##########################################################################

Full_results_output$P_H_ratio <- abs(Full_results_output$P_est.mean / Full_results_output$H_est.mean)

####################################################
##Run classifier to establish class of each region##
####################################################
##Set qvalue cut-off
critical_value <- 0.05

##Run classifier
Full_results_output$Regulatory_class <- "Conserved/Ambiguous"

for (i in 1:nrow(Full_results_output)) {

if (Full_results_output$P_qvalue[i] > critical_value & Full_results_output$H_qvalue[i] > critical_value & Full_results_output$P_H_qvalue[i] > critical_value){

	Full_results_output$Regulatory_class[i] <- "Conserved/Ambiguous"

} else if (Full_results_output$P_qvalue[i] < critical_value & Full_results_output$H_qvalue[i] < critical_value & Full_results_output$P_H_qvalue[i] > critical_value){

	Full_results_output$Regulatory_class[i] <- "Cis"

} else if (Full_results_output$P_qvalue[i] < critical_value & Full_results_output$H_qvalue[i] > critical_value & Full_results_output$P_H_qvalue[i] < critical_value){

	Full_results_output$Regulatory_class[i] <- "Trans"

} else if (Full_results_output$P_qvalue[i] < critical_value & Full_results_output$H_qvalue[i] < critical_value & Full_results_output$P_H_qvalue[i] < critical_value & Full_results_output$Direction_parent[i] == Full_results_output$Direction_hybrid[i] & Full_results_output$P_H_ratio[i] > 1 & Full_results_output$Direction_hybrid[i] != "Ambig" & Full_results_output$Direction_parent[i] != "Ambig"){

	Full_results_output$Regulatory_class[i] <- "Cis_+_Trans,opposing"

} else if (Full_results_output$P_qvalue[i] < critical_value & Full_results_output$H_qvalue[i] < critical_value & Full_results_output$P_H_qvalue[i] < critical_value & Full_results_output$Direction_parent[i] == Full_results_output$Direction_hybrid[i] & Full_results_output$P_H_ratio[i] < 1 & Full_results_output$Direction_hybrid[i] != "Ambig" & Full_results_output$Direction_parent[i] != "Ambig"){

  Full_results_output$Regulatory_class[i] <- "Cis_+_Trans,same"

} else if (Full_results_output$P_qvalue[i] < critical_value & Full_results_output$H_qvalue[i] < critical_value & Full_results_output$P_H_qvalue[i] < critical_value & Full_results_output$Direction_parent[i] != Full_results_output$Direction_hybrid[i] & Full_results_output$Direction_hybrid[i] != "Ambig" & Full_results_output$Direction_parent[i] != "Ambig"){

	Full_results_output$Regulatory_class[i] <- "Cis_*_Trans"

} else if (Full_results_output$P_qvalue[i] > critical_value & Full_results_output$H_qvalue[i] < critical_value & Full_results_output$P_H_qvalue[i] < critical_value){

	Full_results_output$Regulatory_class[i] <- "Compensatory"
}
}

##Run classifier for opposing and same
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

nrow(subset(Full_results_output, Full_results_output$Direction == "Opposing"))
nrow(subset(Full_results_output, Full_results_output$Direction == "Reinforcing"))


#### Compute % CIS and TRANS ####
Full_results_output$trans_reg_diff <- Full_results_output$P_est.mean - Full_results_output$H_est.mean
Full_results_output$perc_cis <- (abs(Full_results_output$H_est.mean)/(abs(Full_results_output$H_est.mean) + abs(Full_results_output$trans_reg_diff))) * 100


##################################
##Write out full results to file##
##################################head

write.table(Full_results_output, file = "~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_intergenic.txt", sep = "\t", row.names = F, quote = F)

##########################
##Generate summary plots##
##########################

perc_cis <- ggplot(Full_results_output, aes(perc_cis)) +
geom_density(color="darkblue", fill="lightblue") +
theme_main() +
xlab("Percent Cis") +
ggtitle("D.mel,ZHR - D.mel,Z30 INTERGENIC Chromatin accessiblity divergence")
ggsave(perc_cis, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/perc_cis_ZHR_Z30_sub_ATAC_CPM_20min_intergenic.pdf", width = 15, height = 15)


cis_trans_ATAC_CPM <- Full_results_output %>%

ggplot(aes(x = P_est.mean, y = H_est.mean)) +
geom_point(alpha = 0.3) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + theme_main() + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-2, 2) + ylim(-2, 2)
ggsave(cis_trans_ATAC_CPM, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/cis_trans_ZHR_Z30_ATAC_sub_CPM_20min_INTERGENIC.pdf", width = 15, height = 15)


############ DENSITY PLOTS  ######
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


Full_results_output$density <- NA

Full_results_output$density[Full_results_output$Regulatory_class == "Cis"] <- get_density(Full_results_output$H_est.mean[Full_results_output$Regulatory_class == "Cis"], Full_results_output$P_est.mean[Full_results_output$Regulatory_class == "Cis"], n = 500)

Full_results_output$density[Full_results_output$Regulatory_class == "Trans"] <- get_density(Full_results_output$H_est.mean[Full_results_output$Regulatory_class == "Trans"], Full_results_output$P_est.mean[Full_results_output$Regulatory_class == "Trans"], n = 500)

Full_results_output$density[Full_results_output$Regulatory_class == "Cis_*_Trans"] <- get_density(Full_results_output$H_est.mean[Full_results_output$Regulatory_class == "Cis_*_Trans"], Full_results_output$P_est.mean[Full_results_output$Regulatory_class == "Cis_*_Trans"], n = 500)

Full_results_output$density[Full_results_output$Regulatory_class == "Cis_+_Trans,same"] <- get_density(Full_results_output$H_est.mean[Full_results_output$Regulatory_class == "Cis_+_Trans,same"], Full_results_output$P_est.mean[Full_results_output$Regulatory_class == "Cis_+_Trans,same"], n = 500)

Full_results_output$density[Full_results_output$Regulatory_class == "Cis_+_Trans,opposing"] <- get_density(Full_results_output$H_est.mean[Full_results_output$Regulatory_class == "Cis_+_Trans,opposing"], Full_results_output$P_est.mean[Full_results_output$Regulatory_class == "Cis_+_Trans,opposing"], n = 500)

Full_results_output$density[Full_results_output$Regulatory_class == "Compensatory"] <- get_density(Full_results_output$H_est.mean[Full_results_output$Regulatory_class == "Compensatory"], Full_results_output$P_est.mean[Full_results_output$Regulatory_class == "Compensatory"], n = 500)

Full_results_output$density[Full_results_output$Regulatory_class == "Conserved/Ambiguous"] <- get_density(Full_results_output$H_est.mean[Full_results_output$Regulatory_class == "Conserved/Ambiguous"], Full_results_output$P_est.mean[Full_results_output$Regulatory_class == "Conserved/Ambiguous"], n = 500)



A <- Full_results_output %>% subset(Regulatory_class == "Cis") %>%

ggplot(aes(x = P_est.mean, y = H_est.mean)) + geom_point(alpha = 0.9) +
geom_point(aes(P_est.mean, H_est.mean, col = density)) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + theme_main() + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-3, 3) + ylim(-3, 3) +  guides(col = F) + ggtitle("Cis")

B <- Full_results_output %>% subset(Regulatory_class == "Trans") %>%

ggplot(aes(x = P_est.mean, y = H_est.mean)) + geom_point(alpha = 0.9) +
geom_point(aes(P_est.mean, H_est.mean, col = density)) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + theme_main() + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-3, 3) + ylim(-3, 3) +  guides(col = F) + ggtitle("Trans")

C <- Full_results_output %>% subset(Regulatory_class == "Cis_*_Trans") %>%

ggplot(aes(x = P_est.mean, y = H_est.mean)) + geom_point(alpha = 0.9) +
geom_point(aes(P_est.mean, H_est.mean, col = density)) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + theme_main() + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-3, 3) + ylim(-3, 3) +  guides(col = F) + ggtitle("Cis * Trans")

D <- Full_results_output %>% subset(Regulatory_class == "Cis_+_Trans,same") %>%

ggplot(aes(x = P_est.mean, y = H_est.mean)) + geom_point(alpha = 0.9) +
geom_point(aes(P_est.mean, H_est.mean, col = density)) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + theme_main() + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-3, 3) + ylim(-3, 3) +  guides(col = F) + ggtitle("Cis + Trans,same")


E <- Full_results_output %>% subset(Regulatory_class == "Cis_+_Trans,opposing") %>%

ggplot(aes(x = P_est.mean, y = H_est.mean)) + geom_point(alpha = 0.9) +
geom_point(aes(P_est.mean, H_est.mean, col = density)) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + theme_main() + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-3, 3) + ylim(-3, 3) +  guides(col = F) + ggtitle("Cis + Trans,opposing")

G <- Full_results_output %>% subset(Regulatory_class == "Compensatory") %>%

ggplot(aes(x = P_est.mean, y = H_est.mean)) + geom_point(alpha = 0.9) +
geom_point(aes(P_est.mean, H_est.mean, col = density)) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + theme_main() + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-3, 3) + ylim(-3, 3) +  guides(col = F) + ggtitle("Compensatory")

H <- Full_results_output %>% subset(Regulatory_class == "Conserved/Ambiguous") %>%

ggplot(aes(x = P_est.mean, y = H_est.mean)) + geom_point(alpha = 0.9) +
geom_point(aes(P_est.mean, H_est.mean, col = density)) +
geom_hline(yintercept = 0, linetype = "dashed") +
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + theme_main() + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-3, 3) + ylim(-3, 3) +  guides(col = F) + ggtitle("Conserved/Ambiguous")

facet_all <- plot_grid(A, B, C, D, E, G, H)

ggsave(facet_all, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/cis_trans_ZHR_Z30_ATAC_CPM_20min_ALL_reg_classes_facet_INTERGENIC.pdf", width = 15, height = 15)

facet_cis_trans <- plot_grid(C, D, E, G)

ggsave(facet_cis_trans, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/cis_trans_ZHR_Z30_ATAC_CPM_20min_cis_trans_reg_classes_facet_INTERGENIC.pdf", width = 15, height = 15)
