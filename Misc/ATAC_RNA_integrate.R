# generate bedtools closest file by finding closest gene coordinates (or first exon) OF EXPRESSED GENES IN RNA-SEQ
## generated input files by extracting coordinates from RNA and ATAC files (see below) and sorting by coordinate

# attach first exon coords to RNA data
RNA <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_100max.txt", header = T)
first_exon_coords <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_gene_coords_first_exon.bed", header = F)
colnames(first_exon_coords) <- c("chrom_exon1", "start_exon1", "end_exon1", "Gene")
RNA_coords <- join_all(list(RNA, first_exon_coords), by = "Gene", type = "left")
write.table(RNA_coords, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_100max_exon1_coords.txt", row.names = F, quote = F, sep = "\t")

# rearrange RNA_coords in unix to just be coordinates and in BED format for below

bedtools closest \
-a /Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/ZHR_Z30_Full_results_output_ALL_classes_ONLY_coords_sorted.bed \
-b /Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_100max_coords_only_sorted.bed \
| uniq \
> /Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_RNA_comp/ZHR_Z30_ATAC_RNA_closest.bed

# read in data from bedtools
ATAC <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/ZHR_Z30_Full_results_output_ALL_classes.txt", header = T)
RNA <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_100max.txt", header = T)
ATAC_RNA_closest  <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_RNA_comp/ZHR_Z30_ATAC_RNA_closest.bed", header = F, sep = "\t")
## clean up ATAC and RNA to only keep necessary columns
# prepare ATAC and RNA files to merge
RNA_minimal <- RNA[,c(1,28:40,43)]
ATAC_minimal <- ATAC[,c(1:4,19:31,34,35)]

## select classes
ATAC_minimal_intra_inter <-  ATAC_minimal[ATAC_minimal$class == "inter" | ATAC_minimal$class == "intra",]
ATAC_minimal_start_end <- ATAC_minimal[ATAC_minimal$class == "start" | ATAC_minimal$class == "end",]


colnames(RNA_minimal) <- paste(colnames(RNA_minimal), "RNA", sep = "_")
colnames(RNA_minimal)[1] <- "Gene"

## clean up closest file to only keep needed columnns and create a locus key
ATAC_RNA_closest$Paste_locus <- paste(ATAC_RNA_closest$V1, ATAC_RNA_closest$V2, ATAC_RNA_closest$V3, sep = "_")
ATAC_RNA_closest <- ATAC_RNA_closest[,c(5,7,8)]
colnames(ATAC_RNA_closest) <- c("beg_first_exon", "Gene", "Paste_locus")

# join closest gene and exon coordinate to ATAC_inter_intra file
ATAC_minimal_intra_inter_closest_gene <- join_all(list(ATAC_minimal_intra_inter, ATAC_RNA_closest), by = "Paste_locus", type = "full") %>% unique() %>% na.omit()

# join gene expression information
RNA_minimal_closest_locus <- join_all(list(RNA_minimal, ATAC_minimal_intra_inter_closest_gene), by = "Gene", type = "full") %>% unique() %>% na.omit()



# prepare txSTart and End files
start_end <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_all_uniq", header = T)
gene_conversions <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/Dmel_geneID_conversion_table.txt", header = F)
gene_conversions <- gene_conversions[,c(1,5)]
colnames(gene_conversions) <- c("Gene", "ID")

start_minimal <- start_end[,c(1,3,7)]
start_minimal$start_up500 <- start_minimal[,2] - 500
start_minimal$start_down500 <- start_minimal[,2] + 500
start_minimal$Paste_locus <- paste(start_minimal[,1], start_minimal[,4], start_minimal[,5], sep = "_")
start_minimal <- start_minimal[,c(6,3)]
colnames(start_minimal) <- c("Paste_locus", "Gene")

start_key <- join_all(list(start_minimal, gene_conversions), type = "left", by = "Gene")


end_minimal <- start_end[,c(1,4,7)]
end_minimal$end_up500 <- end_minimal[,2] - 500
end_minimal$end_down500 <- end_minimal[,2] + 500
end_minimal$Paste_locus <- paste(end_minimal[,1], end_minimal[,4], end_minimal[,5], sep = "_")
end_minimal <- end_minimal[,c(6,3)]
colnames(end_minimal) <- c("Paste_locus", "Gene")

end_key <- join_all(list(end_minimal, gene_conversions), type = "left", by = "Gene")

start_end_key <- rbind(start_key, end_key)

# join associated gene with paste_locus from start and end files

ATAC_minimal_start_end_gene <- join_all(list(ATAC_minimal_start_end, start_end_key), type = "left", by = "Paste_locus") %>% unique() %>% na.omit()

colnames(ATAC_minimal_start_end_gene)[20] <- "Name"
colnames(ATAC_minimal_start_end_gene)[21] <- "Gene"

ATAC_minimal_start_end_gene_RNA_int <- join_all(list(RNA_minimal, ATAC_minimal_start_end_gene), type = "full", by = "Gene") %>% unique() %>% na.omit()
ATAC_minimal_start_end_gene_RNA_int$Name <- NULL

# calculate basic metrics with integrated matrix
## distance to first exon
RNA_minimal_closest_locus$distance_to_exon1 <- abs(RNA_minimal_closest_locus[,17] - RNA_minimal_closest_locus[,36])
write.table(RNA_minimal_closest_locus, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_RNA_comp/ZHR_Z30_ATAC_RNA_integrated_minimal.txt", row.names = F, quote = F, sep = "\t")


RNA_minimal_closest_locus$distance_bin <- "7"

for (i in 1:nrow(RNA_minimal_closest_locus)) {

if (RNA_minimal_closest_locus$distance_to_exon1[i] < 5000){

	RNA_minimal_closest_locus$distance_bin[i] <- "1"

} else if (RNA_minimal_closest_locus$distance_to_exon1[i] > 5000 & RNA_minimal_closest_locus$distance_to_exon1[i] < 10000){

	RNA_minimal_closest_locus$distance_bin[i] <- "2"

} else if (RNA_minimal_closest_locus$distance_to_exon1[i] > 10000 & RNA_minimal_closest_locus$distance_to_exon1[i] < 15000){

	RNA_minimal_closest_locus$distance_bin[i] <- "3"

} else if (RNA_minimal_closest_locus$distance_to_exon1[i] > 15000 & RNA_minimal_closest_locus$distance_to_exon1[i] < 20000){

	RNA_minimal_closest_locus$distance_bin[i] <- "4"

} else if (RNA_minimal_closest_locus$distance_to_exon1[i] > 20000 & RNA_minimal_closest_locus$distance_to_exon1[i] < 25000){

	RNA_minimal_closest_locus$distance_bin[i] <- "5"

} else if (RNA_minimal_closest_locus$distance_to_exon1[i] > 25000 & RNA_minimal_closest_locus$distance_to_exon1[i] < 30000){

RNA_minimal_closest_locus$distance_bin[i] <- "6"

}
}


##### assign txEnd and txStart their genes
