## correlation heatmap

# ATAC
ALL <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes.txt", header = T)

ALL_data <- ALL[,c(5:16)]

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }


  reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  }

cormat <- round(cor(ALL_data, method = "pearson"),2)
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()
 ggheatmap +
geom_text(aes(Var2, Var1, label = value), color = "grey", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.6, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))


# RNA
RNA <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_1000max.txt", header = T)

ALL_data <- RNA[,c(2:25)]

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }


  reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  }

cormat <- round(cor(ALL_data, method = "pearson"),2)
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()
 ggheatmap +
geom_text(aes(Var2, Var1, label = value), color = "grey", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.6, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))


## variation!!

ZHR_Z30_regions_ALL_SNPs <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_SNPs_intersectALL.bed", header = F)
colnames(ZHR_Z30_regions_ALL_SNPs) <- c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2",
"pos_start", "pos_end", "P_est.mean", "P_p_value","H_est.mean","H_p_value","H_P_est","H_P_p_value","P_qvalue","H_qvalue","P_H_qvalue","Direction_parent","Direction_hybrid",
"P_H_ratio","Regulatory_class","trans_reg_diff","perc_cis","class", "chrom_snp", "start_snp","end_snp","overlap_snp")

## clean up and add column with number of total SNPs in region
ZHR_Z30_regions_ALL_SNPs <- ZHR_Z30_regions_ALL_SNPs[,c(1:(ncol(ZHR_Z30_regions_ALL_SNPs) - 4))]

ZHR_Z30_regions_ALL_SNPs <- ddply(ZHR_Z30_regions_ALL_SNPs,.(chrom, start, end, Paste_locus, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, HYB_1_P1, HYB_2_P1, HYB_3_P1, HYB_1_P2, HYB_2_P2, HYB_3_P2,
pos_start, pos_end, P_est.mean, P_p_value,H_est.mean,H_p_value,H_P_est,H_P_p_value,P_qvalue,H_qvalue,P_H_qvalue,Direction_parent,Direction_hybrid,
P_H_ratio,Regulatory_class,trans_reg_diff,perc_cis,class),nrow)

ZHR_Z30_regions_ALL_SNPs_final <- ZHR_Z30_regions_ALL_SNPs[,c(4,ncol(ZHR_Z30_regions_ALL_SNPs))]
colnames(ZHR_Z30_regions_ALL_SNPs_final) <- c("Paste_locus", "num_SNPs_total")

ZHR_Z30_regions_ALL_SNPs_final$snp_perc <- (ZHR_Z30_regions_ALL_SNPs_final$num_SNPs_total/1000) * 100

## get CDS snps
bedtools intersect -wao -a zhr_z30_SNPs_dm6_genome.bed -b dm6_CDS_lengths.bed | uniq > zhr_z30_SNPs_CDS_lengths.bed

CDS_snps <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/zhr_z30_SNPs_CDS_lengths.bed", header = F)
CDS_snps <- CDS_snps[CDS_snps$V8 != 0,]
CDS_snps <- CDS_snps[,c(4:7)]
colnames(CDS_snps) <- c("chrom", "start", "end", "length")
CDS_snps <- ddply(CDS_snps,.(chrom, start, end, length),nrow)
colnames(CDS_snps)[5] <- "snp_num"
CDS_snps$length <- CDS_snps$length %>% as.numeric()
CDS_snps$snp_num <- CDS_snps$snp_num %>% as.numeric()
CDS_snps$snp_perc <- (CDS_snps$snp_num/CDS_snps$length)*100


ggplot() +
geom_density(data=CDS_snps,aes(x=snp_perc), fill = "purple", alpha = 0.5) +
geom_density(data=ZHR_Z30_regions_ALL_SNPs_final,aes(x=snp_perc), fill = "green", alpha = 0.5) +
xlim(0,25) +
xlab("Percent sequence divergence")

# median sequence divergence
## CDS = 0.79%
## peaks = 1%

#### relationship between sequence variation and accessibility variation ###
ZHR_Z30_regions_ALL_SNPs <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/ATAC_seq_datafiles/Bayes_test_outputs/Full_results_output_ZHR_Z30_ATAC_20min_centered1000_classes_SNPs_intersectALL.bed", header = F)

colnames(ZHR_Z30_regions_ALL_SNPs) <- c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2",
"pos_start", "pos_end", "P_est.mean", "P_p_value","H_est.mean","H_p_value","H_P_est","H_P_p_value","P_qvalue","H_qvalue","P_H_qvalue","Direction_parent","Direction_hybrid",
"P_H_ratio","Regulatory_class","trans_reg_diff","perc_cis","class", "chrom_snp", "start_snp","end_snp","overlap_snp")

ZHR_Z30_regions_ALL_SNPs <- ZHR_Z30_regions_ALL_SNPs[,c(1:(ncol(ZHR_Z30_regions_ALL_SNPs) - 4))]

ZHR_Z30_regions_ALL_SNPs <- ddply(ZHR_Z30_regions_ALL_SNPs,.(chrom, start, end, Paste_locus, P1_1, P1_2, P1_3, P2_1, P2_2, P2_3, HYB_1_P1, HYB_2_P1, HYB_3_P1, HYB_1_P2, HYB_2_P2, HYB_3_P2,
pos_start, pos_end, P_est.mean, P_p_value,H_est.mean,H_p_value,H_P_est,H_P_p_value,P_qvalue,H_qvalue,P_H_qvalue,Direction_parent,Direction_hybrid,
P_H_ratio,Regulatory_class,trans_reg_diff,perc_cis,class),nrow)

colnames(ZHR_Z30_regions_ALL_SNPs)[ncol(ZHR_Z30_regions_ALL_SNPs)] <- "snp_num"

ggplot(ZHR_Z30_regions_ALL_SNPs) +
geom_point(aes(x=abs(H_est.mean), y=snp_num))

test <- ZHR_Z30_regions_ALL_SNPs

test$hyb_cat <- NA
for (i in 1:nrow(test)) {

if (abs(test$H_est.mean[i]) < 0.128951){

	test$hyb_cat[i] <- "0.128951"

} else if (abs(test$H_est.mean[i]) > 0.128951 & abs(test$H_est.mean[i]) < 0.247514){

	test$hyb_cat[i] <- "0.128951_0.247514"

} else if (abs(test$H_est.mean[i]) > 0.247514 & abs(test$H_est.mean[i]) < 0.548344){

	test$hyb_cat[i] <- "0.247514_0.548344"

} else if (abs(test$H_est.mean[i]) > 0.548344){

	test$hyb_cat[i] <- "0.548344_above"

}
}

ggplot(test, aes(x=abs(H_est.mean)))+
geom_density(fill = "purple", alpha = 0.5)+
theme_main()+
xlab("Hybrid allelic Divergence")+
geom_vline(xintercept=c(0.128951,0.247514,0.548344), linetype = "dotted")


ggplot(test, aes(x=hyb_cat, y=snp_num, color=hyb_cat)) +
geom_boxplot(notch=T, fill="purple") +
theme_main() +
ylab("SNP number") +
xlab("Hybrid div") +
scale_color_discrete(guide=F)

test$class <- factor(test$class,levels = c("start", "end", "inter", "intra"))

ggplot(test, aes(x=class, y=snp_num, fill=class)) +
geom_boxplot(notch=T) +
theme_main() +
ylab("SNP number") +
xlab("Hybrid div") +
scale_color_discrete(guide=F)
