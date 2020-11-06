######### Analyses on Harrison lab dataset for Grh ###########

# QUESTION: is there a correlation between chromatin accessibility and Grh binding intensity?
## Going to take NarrowPeak ATAC-seq file with peak ranges and scores, and calculate the number of reads overlapping each range from Grh ChIP dataset

# go to wd with bam and bed file
cd /Users/henryertl/Documents/Wittkopp_lab/pop_genomics_reg_regions/Grh/embryo_harrison_paper_data

# ensure bam is indexed
samtools index G2-3hr_Grh_ChIP_dm6.bam

# run bedtools multicov
bedtools multicov -bams 2-3hr_Grh_ChIP_dm6.bam -bed ATAC_stage5_peaks.bed > 2_3HR_Grh_ATAC.bed

## load into R and plot
######### R script ##########
#libraries
library(ggplot2);

# read in datasets & add columns names
Grh_ChIP_2_3hrs <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/pop_genomics_reg_regions/Grh/embryo_harrison_paper_data/2_3HR_Grh_ATAC.bed", header = F)
colnames(Grh_ChIP_2_3hrs) <- c("chrom", "start", "end", "name", "ATAC_score", "Grh_reads")
Grh_ChIP_2_3hrs <- as.data.frame(Grh_ChIP_2_3hrs);

# plot
ggplot(Grh_ChIP_2_3hrs) +
  geom_point(aes(x=ATAC_score, y=Grh_reads))


bedtools intersect -wa -a 2_3HR_Grh_ATAC.bed -b 2-3hr_GRH_peaks.bed > 2_3HR_Grh_ATAC_only_Grh_peaks.bed
