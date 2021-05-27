## goal: take output from Define_universal_peakset.sh and categorize into and intergenic, intragenic & create txStart and txEnd BED files

### downloaded txStart and txEnd coords from UCSC browswer and sorted uniq

# TXSTART AND END -- USE EXISTING DF FROM UCSC TO DEFINE +/- 500BP TXSTART AND END

#setwd to folder with cloned 'Integrative_AS_genomics' repository
setwd("/Users/wittkopp_member/Code/")

# read in dataframe and add locus ID
df <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txStart_End_coords.txt", header = T)
colnames(df)[1] <- "chrom"
df$paste_locus <- paste(df$chrom, df$txStart, df$txEnd, sep = "_")

# add txStart +- 500bp and add txEnd +- 500bp
df$txStart_plus500 <- df$txStart + 500
df$txStart_minus500 <- df$txStart - 500
df$txStart_plus1 <- df$txStart + 1

df$txEnd_plus500 <- df$txEnd + 500
df$txEnd_minus500 <- df$txEnd - 500
df$txEnd_plus1 <- df$txEnd + 1


# make final dataframes and write tables
df_txStart <- df[,c("chrom", "txStart_minus500", "txStart_plus500", "paste_locus")]
df_txEnd <- df[,c("chrom", "txEnd_minus500", "txEnd_plus500", "paste_locus")]
df_txStart_End <- df[,c("chrom", "txStart", "txEnd", "paste_locus")]
df_CDS <- df[,c("chrom", "cdsStart", "cdsEnd", "paste_locus")]
df_txStart_only <- df[,c("chrom", "txStart", "txStart_plus1", "paste_locus")]
df_txEnd_only <- df[,c("chrom", "txEnd", "txEnd_plus1", "paste_locus")]

write.table(df_txStart, "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txStart.txt", sep="\t", quote = F, row.names=F)
write.table(df_txStart_End, "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txStart_End.txt", sep="\t", quote = F, row.names=F)
write.table(df_txStart_only, "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txStart_only.txt", sep="\t", quote = F, row.names=F)
write.table(df_txEnd, "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txEnd.txt", sep="\t", quote = F, row.names=F)
write.table(df_txEnd_only, "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txEnd_only.txt", sep="\t", quote = F, row.names=F)
write.table(df_CDS, "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_CDS.txt", sep="\t", quote = F, row.names=F)
write.table(df, "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_all_uniq_processed.txt", sep="\t", quote = F, row.names=F)

# INTERGENIC AND INTRAGENIC CLASSES -- REMOVE ALL PEAKS OVERLAPPING WITH TXSTART AND END

####  BEDTOOLS INTERSECT COMMANDS #####
## Intergenic -- non-overlap with -500 txStart and +500 txEnd
bedtools intersect -v -wa \
-a ZHR_TSIM_ATAC_20min1000max.bed \
-b /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txStart_End.bed \
| uniq \
> /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_TSIM_ATAC_20min1000max_intergenic_NO_TxStart_End.bed

## Intragenic -- overlap with -500 txStart and +500 txEnd
bedtools intersect -wa \
-a ZHR_TSIM_ATAC_20min1000max.bed \
-b /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txStart_End.bed \
| uniq \
> /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/temp.bed

## BUT non-overlap with +/- 500 txStart nor +/- 500 txEnd
bedtools intersect -v -wa \
-a /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/temp.bed \
-b /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txStart.bed \
| uniq \
> /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/temp2.bed

bedtools intersect -v -wa \
-a /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/temp2.bed \
-b /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_txEnd.bed \
| uniq \
> /Users/henryertl/Documents/Devs/Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_TSIM_ATAC_20min1000max_intragenic_NO_TxStart_End.bed


#### Make centered tables from intra- and inter-genic files
# read in inter and intra BEDs from BEDtools to R
inter <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_TSIM_analysis_regions_intergenic_NO_TxStart_End.bed", sep = "\t", header = F)
intra <- read.delim("/./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_TSIM_analysis_regions_intragenic_NO_TxStart_End.bed", sep = "\t", header = F)

# center +/- 500 to compare with txStart and txEnd
inter$length <- inter$V3 - inter$V2
inter_new <- inter[,c(1, 2, 3, 4)]
inter_new$V2 <- (ceiling(inter_new$V2 + (inter_new$length/2))) - 500
inter_new$V3 <- (ceiling(inter_new$V2 + (inter_new$length/2))) + 500
colnames(inter_new) <- c("chrom", "start", "end", "length")
inter_new$length <- NULL
write.table(inter_new, "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_TSIM_analysis_regions_intergenic_NO_TxStart_End_center1000.txt", sep="\t", quote = F, row.names=F)

intra$length <- intra$V3 - intra$V2
intra_new <- intra[,c(1, 2, 3, 4)]
intra_new$V2 <- (ceiling(intra_new$V2 + (intra_new$length/2))) - 500
intra_new$V3 <- (ceiling(intra_new$V2 + (intra_new$length/2))) + 500
colnames(intra_new) <- c("chrom", "start", "end", "length")
intra_new$length <- NULL
write.table(intra_new, "./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/ZHR_Z30_TSIM_analysis_regions_intragenic_NO_TxStart_End_center1000.txt", sep="\t", quote = F, row.names=F)
