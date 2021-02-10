# load libraries and functions
library(stringr)
library(dplyr)
library(plyr)
library(ggplot2)

theme_main <- function() {
  theme_bw() +
  theme(
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  axis.text = element_text(size = 40),
  axis.title = element_text(size = 40),
  strip.text = element_text(size = 30),
  legend.text= element_text(size = 40),
  legend.title = element_text(size = 30),
  plot.title = element_text(size = 25, face = "bold")
)
}

#######################################################################
###### Generate base plots - % cis accumulation and differences #######
#######################################################################

# Read in Bayes output files
C_R <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/Cis_trans_bayes_output_C_R_RNA_CPM.txt", header = T)
C_P <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/Cis_trans_bayes_output_C_P_RNA_CPM.txt", header = T)
C_M <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/Cis_trans_bayes_output_C_M_RNA_CPM.txt", header = T)
C_B <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/Cis_trans_bayes_output_C_B_RNA_CPM.txt", header = T)

# Generate main % cis dataframe
## Assign IDs to properly plot
C_R$ID <- '1'
C_P$ID <- '2'
C_M$ID <- '3'
C_B$ID <- '4'

## Remove genes with conserved/ambiguous expression divergence
C_R_no_cons <- C_R[C_R$Regulatory_class != "Conserved/Ambiguous",]
C_P_no_cons <- C_P[C_P$Regulatory_class != "Conserved/Ambiguous",]
C_M_no_cons <- C_M[C_M$Regulatory_class != "Conserved/Ambiguous",]
C_B_no_cons <- C_B[C_B$Regulatory_class != "Conserved/Ambiguous",]

## Reformat
C_R_perc_cis <- C_R_no_cons[, c("ID", "perc_cis")]
C_P_perc_cis <- C_P_no_cons[, c("ID", "perc_cis")]
C_M_perc_cis <- C_M_no_cons[, c("ID", "perc_cis")]
C_B_perc_cis <- C_B_no_cons[, c("ID", "perc_cis")]
Full_results_ALL_yeast <- rbind(C_R_perc_cis, C_P_perc_cis, C_B_perc_cis, C_M_perc_cis) %>% as.data.frame()

# Plot % cis across species
ALL_perc_cis <- ggplot(Full_results_ALL_yeast, aes(x=ID, y=perc_cis, fill=ID)) +
geom_boxplot(notch=TRUE) +
scale_x_discrete(labels=c("Sc-Sc","Sc-Sp","Sc-Sm","Sc-Sb")) +
xlab("") +
ylab("Percent cis") +
scale_fill_discrete(guide=FALSE) +
theme_main()
  ggsave(ALL_perc_cis, file="/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/evx035f3p_Supp/Figures/Perc_cis_ALL.pdf", width = 12, height = 15)

# Plot distribution of ∆ cis by divergence for C_M and C_B
## Merge two dataframes
C_B_M <- merge(C_M, C_B, by='Gene')

## Calculate ∆ cis
C_B_M$perc_cis_diff <- log2(C_B_M$perc_cis.y/C_B_M$perc_cis.x)

## Remove genes with conserved/ambiguous expression divergence
C_B_M <- C_B_M[C_B_M$Regulatory_class.x != "Conserved/Ambiguous" & C_B_M$Regulatory_class.y != "Conserved/Ambiguous" ,]

# Plot ∆ cis
ggplot(C_B_M, aes(x=perc_cis_diff)) +
geom_density(fill = "#E69F00") +
theme_main() +
xlab("∆ cis by divergence Sc/b-Sc/m") +
geom_vline(xintercept=0)
  ggsave(S, file = "/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast//evx035f3p_Supp/Figures/perc_cis_divergence_distribution_Sc_b_m_RNA.pdf", width = 10, height = 15)


Full_results_comb_coords <- C_B_M

# Assign percentile expression
## get summary stats on expression for P1_1.x
stat_sum <- summary(Full_results_comb_coords$P1_1.x)

## Assign expression quartile to each gene
Full_results_comb_coords$Expression_quantile <- "NA"

for (i in 1:nrow(Full_results_comb_coords)) {

if (Full_results_comb_coords$P1_1.x[i] <= as.numeric(stat_sum[2])){

	Full_results_comb_coords$Expression_quantile[i] <- "1Q"

} else if (Full_results_comb_coords$P1_1.x[i] >= as.numeric(stat_sum[2]) & Full_results_comb_coords$P1_1.x[i] <= as.numeric(stat_sum[3])){

	Full_results_comb_coords$Expression_quantile[i] <- "2Q"

} else if (Full_results_comb_coords$P1_1.x[i] >= as.numeric(stat_sum[3]) & Full_results_comb_coords$P1_1.x[i] <= as.numeric(stat_sum[5])){

	Full_results_comb_coords$Expression_quantile[i] <- "3Q"

} else if (Full_results_comb_coords$P1_1.x[i] >= as.numeric(stat_sum[5]) & Full_results_comb_coords$P1_1.x[i] <= as.numeric(stat_sum[6])){

	Full_results_comb_coords$Expression_quantile[i] <- "4Q"

}
}


## % cis change with divergence: by expression quantile
O <- ggplot(Full_results_comb_coords, aes(x = Expression_quantile, y = perc_cis_diff, fill = Expression_quantile))+
#geom_violin() +
geom_boxplot(notch=TRUE) +
lims(y=c(-3,3)) +
#ylim(-2.5,2.5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
geom_smooth(method="loess", aes(group=1))
ggsave(O, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Figures/perc_cis_div_change_Express_quartile.pdf", width = 10, height = 15)

# Estimated expression divergence by expression quantile
P <- ggplot(Full_results_comb_coords, aes(x=Expression_quantile, y=abs(P_est.mean.y), fill = Expression_quantile))+
#geom_violin() +
geom_boxplot(notch=TRUE) +
ylim(0,1) +
xlab("") +
ylab("Estimated expression divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE)
ggsave(P, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Figures/Expression_div_Express_quartile.pdf", width = 10, height = 15)

# Expression profilr divergence by expression quantile
x <- Full_results_comb_coords[Full_results_comb_coords$Auto_X == "X",]
auto <- Full_results_comb_coords[Full_results_comb_coords$Auto_X == "Autosome",]
P1_cor_X <- x[, c("P1_1.y", "P2_1.y")]
P1_cor_Auto <- auto[, c("P1_1.y", "P2_1.y")]
cor(P1_cor_X, method = "spearman")
cor(P1_cor_Auto, method = "spearman")


Q1 <- Full_results_comb_coords[Full_results_comb_coords$Expression_quantile == "1Q",]
Q2 <- Full_results_comb_coords[Full_results_comb_coords$Expression_quantile == "2Q",]
Q3 <- Full_results_comb_coords[Full_results_comb_coords$Expression_quantile == "3Q",]
Q4 <- Full_results_comb_coords[Full_results_comb_coords$Expression_quantile == "4Q",]

# P1
P1_cor_Q1 <- Q1[, c("P1_1.x", "P2_1.x")]
P1_cor_Q2 <- Q2[, c("P1_1.x", "P2_1.x")]
P1_cor_Q3 <- Q3[, c("P1_1.x", "P2_1.x")]
P1_cor_Q4 <- Q4[, c("P1_1.x", "P2_1.x")]


cor(P1_cor_Q1, method = "spearman")
cor(P1_cor_Q2, method = "spearman")
cor(P1_cor_Q3, method = "spearman")
cor(P1_cor_Q4, method = "spearman")

# P2
P2_cor_Q1 <- Q1[, c("P1_2.x", "P2_2.x")]
P2_cor_Q2 <- Q2[, c("P1_2.x", "P2_2.x")]
P2_cor_Q3 <- Q3[, c("P1_2.x", "P2_2.x")]
P2_cor_Q4 <- Q4[, c("P1_2.x", "P2_2.x")]

cor(P2_cor_Q1, method = "spearman")
cor(P2_cor_Q2, method = "spearman")
cor(P2_cor_Q3, method = "spearman")
cor(P2_cor_Q4, method = "spearman")

# P3
P3_cor_Q1 <- Q1[, c("P1_3.y", "P2_3.y")]
P3_cor_Q2 <- Q2[, c("P1_3.y", "P2_3.y")]
P3_cor_Q3 <- Q3[, c("P1_3.y", "P2_3.y")]
P3_cor_Q4 <- Q4[, c("P1_3.y", "P2_3.y")]

cor(P3_cor_Q1, method = "spearman")
cor(P3_cor_Q2, method = "spearman")
cor(P3_cor_Q3, method = "spearman")
cor(P3_cor_Q4, method = "spearman")


# Is the decreasing pattern of expression divergence by expression quantile an artifact of lower numbers?
## Look at for just 3Q and 4Q
Full_results_comb_coords_3Q <- subset(Full_results_comb_coords, Full_results_comb_coords$Expression_quantile == "3Q")
Full_results_comb_coords_4Q <- subset(Full_results_comb_coords, Full_results_comb_coords$Expression_quantile == "4Q")
Full_results_comb_coords_3Q_4Q <- rbind(Full_results_comb_coords_3Q, Full_results_comb_coords_4Q)

summary(Full_results_comb_coords_3Q_4Q$P1_1.x)

Full_results_comb_coords_3Q_4Q$High_expression_quantile <- "1Q"

for (i in 1:nrow(Full_results_comb_coords_3Q_4Q)) {

if (Full_results_comb_coords_3Q_4Q$P1_1.x[i] < 395){

	Full_results_comb_coords_3Q_4Q$High_expression_quantile[i] <- "1Q"

} else if (Full_results_comb_coords_3Q_4Q$P1_1.x[i] > 395 & Full_results_comb_coords_3Q_4Q$P1_1.x[i] < 491){

	Full_results_comb_coords_3Q_4Q$High_expression_quantile[i] <- "2Q"

} else if (Full_results_comb_coords_3Q_4Q$P1_1.x[i] > 491 & Full_results_comb_coords_3Q_4Q$P1_1.x[i] < 618){

	Full_results_comb_coords_3Q_4Q$High_expression_quantile[i] <- "3Q"

} else if (Full_results_comb_coords_3Q_4Q$P1_1.x[i] > 618 & Full_results_comb_coords_3Q_4Q$P1_1.x[i] < 1219){

	Full_results_comb_coords_3Q_4Q$High_expression_quantile[i] <- "4Q"

}
}

# Estimated expression divergence by expression quantile - 3Q and 4Q only
Q <- ggplot(Full_results_comb_coords_3Q_4Q, aes(x=High_expression_quantile, y=abs(P_est.mean.x), fill = High_expression_quantile))+
#geom_violin() +
geom_boxplot(notch=TRUE) +
ylim(0,0.5) +
xlab("Quartile ranges for upper two quartiles of all data") +
ylab("Estimated expression divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE)
ggsave(Q, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Figures/Expression_div_upper2Q_Express_quartile.pdf", width = 10, height = 15)



## Get summary stats on expression for P1_1.x
stat_sum <- summary(abs(Full_results_comb_coords$P_est.mean.y))

## Assign expression quartile to each gene
Full_results_comb_coords$Expression_change_quantile <- "NA"

for (i in 1:nrow(Full_results_comb_coords)) {

if (Full_results_comb_coords$P_est.mean.y[i] <= as.numeric(stat_sum[2])){

	Full_results_comb_coords$Expression_change_quantile[i] <- "1Q"

} else if (Full_results_comb_coords$P_est.mean.y[i] >= as.numeric(stat_sum[2]) & Full_results_comb_coords$P_est.mean.y[i] <= as.numeric(stat_sum[3])){

	Full_results_comb_coords$Expression_change_quantile[i] <- "2Q"

} else if (Full_results_comb_coords$P_est.mean.y[i] >= as.numeric(stat_sum[3]) & Full_results_comb_coords$P_est.mean.y[i] <= as.numeric(stat_sum[5])){

	Full_results_comb_coords$Expression_change_quantile[i] <- "3Q"

} else if (Full_results_comb_coords$P_est.mean.y[i] >= as.numeric(stat_sum[5]) & Full_results_comb_coords$P_est.mean.y[i] <= as.numeric(stat_sum[6])){

	Full_results_comb_coords$Expression_change_quantile[i] <- "4Q"

}
}

ggplot(Full_results_comb_coords, aes(x = Expression_change_quantile, y = perc_cis_diff, fill = Expression_change_quantile))+
#geom_violin() +
geom_boxplot(notch=TRUE) +
ylim(-2.5,2.5) +
xlab("Expression Quantile") +
ylab("% cis change with divergence") +
theme_main() +
scale_fill_discrete(guide=FALSE) +
geom_smooth(method="loess", aes(group=1))
