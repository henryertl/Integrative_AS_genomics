# generate bedtools closest file by finding closest gene coordinates (or first exon) OF EXPRESSED GENES IN RNA-SEQ
## generated input files by extracting coordinates from RNA and ATAC files (see below) and sorting by coordinate



##################################
## generating files input files ##
##################################

setwd("/Users/henryertl/Documents/Devs/Integrative_AS_genomics")

# attach first exon coords to RNA data
RNA <- read.delim("./AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_100max.txt", header = T)
first_exon_coords <- read.delim("./AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_gene_coords_first_exon.bed", header = F)
colnames(first_exon_coords) <- c("chrom_exon1", "start_exon1", "end_exon1", "Gene")
RNA_coords <- join_all(list(RNA, first_exon_coords), by = "Gene", type = "left")
write.table(RNA_coords, file = "./AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_100max_exon1_coords.txt", row.names = F, quote = F, sep = "\t")

# rearrange RNA_coords in unix to just be coordinates and in BED format for below
bedtools closest \
-a /Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/ZHR_Z30_Full_results_output_ALL_classes_ONLY_coords_sorted.bed \
-b /Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_100max_coords_only_sorted.bed \
| uniq \
> /Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_RNA_comp/ZHR_Z30_ATAC_RNA_closest.bed

########################################
## join files for integrated analyses ##
########################################

# read in data from bedtools
ATAC <- read.delim("./AS_ATAC_RNA_2020_10_1/ATAC_seq/ZHR_Z30_Full_results_output_ALL_classes.txt", header = T)
RNA <- read.delim("./AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/Full_results_output_ZHR_Z30_RNA_20min_100max.txt", header = T)
ATAC_RNA_closest  <- read.delim("./AS_ATAC_RNA_2020_10_1/ATAC_RNA_comp/ZHR_Z30_ATAC_RNA_closest.bed", header = F, sep = "\t")

## clean up ATAC data - keep relevant columns and separate classes
ATAC_minimal <- ATAC[,c(1:4,19:31,34,35)]
ATAC_minimal_intra_inter <-  ATAC_minimal[ATAC_minimal$class == "inter" | ATAC_minimal$class == "intra",]
ATAC_minimal_start_end <- ATAC_minimal[ATAC_minimal$class == "start" | ATAC_minimal$class == "end",]

## clean up RNA data  - keep relevant columns and rename columns in prep for join with ATAC
RNA_minimal <- RNA[,c(1,28:40,43)]
colnames(RNA_minimal) <- paste(colnames(RNA_minimal), "RNA", sep = "_")
colnames(RNA_minimal)[1] <- "Gene"

## clean up closest file - add locus key to join with
ATAC_RNA_closest$Paste_locus <- paste(ATAC_RNA_closest$V1, ATAC_RNA_closest$V2, ATAC_RNA_closest$V3, sep = "_")
ATAC_RNA_closest <- ATAC_RNA_closest[,c(5,7,8)]
colnames(ATAC_RNA_closest) <- c("beg_first_exon", "Gene", "Paste_locus")

# SET 1: integrate closest gene info for inter- and intra-genic ATAC regions
## join closest gene and exon coordinate to inter- and intra-genic regions
ATAC_minimal_intra_inter_closest_gene <- join_all(list(ATAC_minimal_intra_inter, ATAC_RNA_closest), by = "Paste_locus", type = "full") %>% unique() %>% na.omit()
## join gene expression information
RNA_minimal_closest_locus <- join_all(list(RNA_minimal, ATAC_minimal_intra_inter_closest_gene), by = "Gene", type = "full") %>% unique() %>% na.omit()
## OUTPUT 4410 regions with closest gene info



# Set 2: integrate start and end data with corresponding gene and expression info
# prepare txSTart and End files - basically just re-writing the locus key
start_end <- read.delim("./AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/dm6_all_uniq", header = T)
gene_conversions <- read.delim("./AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/Dmel_geneID_conversion_table.txt", header = F)
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

#### this is the final key for txStart and txEnd coordinates and gene info
start_end_key <- rbind(start_key, end_key)

# join associated gene with paste_locus from start and end files
ATAC_minimal_start_end_gene <- join_all(list(ATAC_minimal_start_end, start_end_key), type = "left", by = "Paste_locus") %>% unique() %>% na.omit()
colnames(ATAC_minimal_start_end_gene)[20] <- "Name"
colnames(ATAC_minimal_start_end_gene)[21] <- "Gene"

# join gene expression info
ATAC_minimal_start_end_gene_RNA_int <- join_all(list(RNA_minimal, ATAC_minimal_start_end_gene), type = "full", by = "Gene") %>% unique() %>% na.omit()
ATAC_minimal_start_end_gene_RNA_int$Name <- NULL
ATAC_minimal_start_end_gene_RNA_int$beg_first_exon <- 0

## OUTPUT: 3768 txStart/End regions

### join 2 classes ###
all_classes_integrated <- rbind(ATAC_minimal_start_end_gene_RNA_int, RNA_minimal_closest_locus)

### ANALYSES ###

# Categorize

all_classes_integrated$div_category <- "NA"

for (i in 1:nrow(all_classes_integrated)) {

if (all_classes_integrated$P_qvalue_RNA[i] < 0.05 && all_classes_integrated$P_qvalue[i] < 0.05 && sign(all_classes_integrated$P_est.mean_RNA[i]) == sign(all_classes_integrated$P_est.mean[i])){

	all_classes_integrated$div_category[i] <- "AccD_ExpD_S"

} else if (all_classes_integrated$P_qvalue_RNA[i] < 0.05 && all_classes_integrated$P_qvalue[i] < 0.05 && sign(all_classes_integrated$P_est.mean_RNA[i]) != sign(all_classes_integrated$P_est.mean[i])){

	all_classes_integrated$div_category[i] <- "AccD_ExpD_O"

} else if (all_classes_integrated$P_qvalue_RNA[i] > 0.05 && all_classes_integrated$P_qvalue[i] < 0.05){

	all_classes_integrated$div_category[i] <- "AccD_ExpC"

} else if (all_classes_integrated$P_qvalue_RNA[i] < 0.05 && all_classes_integrated$P_qvalue[i] > 0.05){

	all_classes_integrated$div_category[i] <- "AccC_ExpD"

} else if (all_classes_integrated$P_qvalue_RNA[i] > 0.05 && all_classes_integrated$P_qvalue[i] > 0.05){

  all_classes_integrated$div_category[i] <- "AccC_ExpC"

}
}

all_classes_integrated$distance_to_exon1 <- abs(all_classes_integrated[,18] - all_classes_integrated[,35])

write.table(all_classes_integrated, file = "./AS_ATAC_RNA_2020_10_1/ATAC_RNA_comp/ZHR_Z30_ATAC_RNA_integrated_minimal.txt", row.names = F, quote = F, sep = "\t")

all_classes_integrated <- read.delim("./AS_ATAC_RNA_2020_10_1/ATAC_RNA_comp/ZHR_Z30_ATAC_RNA_integrated_minimal.txt", header = T)

O <- all_classes_integrated %>%
ggplot(aes(x=P_est.mean, y=P_est.mean_RNA, color=div_category)) +
geom_point() +
theme_main() +
scale_color_discrete(name = "Divergence category",
labels = c("Acc cons & Exp cons", "Acc cons & Exp div", "Acc div & Exp cons", "Acc div & Exp div opposite", "Acc div & Exp div same ")) +
theme(legend.text= element_text(size = 10),
legend.title = element_text(size = 10)) +
xlab("Estimated accessibility divergence") +
ylab("Estimated expression divergennce") +
ylim(-2.5,2.5) +
xlim(-2.5,2.5) +
ggtitle("Accessibility vs Expression divergence - integrated categories")
	ggsave(O, file = "./Figures/Acc_vs_Exp_categories_ZHR_Z30.pdf")

## without guide
	G <- all_classes_integrated %>%
	ggplot(aes(x=P_est.mean, y=P_est.mean_RNA, color=div_category)) +
	geom_point() +
	theme_main() +
	scale_color_discrete(guide = F) +
	xlab("Estimated accessibility divergence") +
	ylab("Estimated expression divergennce") +
	ylim(-2.5,2.5) +
	xlim(-2.5,2.5) +
	ggtitle("Accessibility vs Expression divergence - integrated categories")

	X <- all_classes_integrated[all_classes_integrated$class == "start",] %>%
	ggplot(aes(x=P_est.mean, y=P_est.mean_RNA, color=div_category)) +
	geom_point() +
	theme_main() +
	scale_color_discrete(guide = F) +
	xlab("Estimated accessibility divergence") +
	ylab("Estimated expression divergennce") +
	ylim(-2.5,2.5) +
	xlim(-2.5,2.5) +
	ggtitle("Accessibility vs Expression divergence - integrated categories")

	Y <- all_classes_integrated[all_classes_integrated$class == "end",] %>%
	ggplot(aes(x=P_est.mean, y=P_est.mean_RNA, color=div_category)) +
	geom_point() +
	theme_main() +
	scale_color_discrete(guide = F) +
	xlab("Estimated accessibility divergence") +
	ylab("Estimated expression divergennce") +
	ylim(-2.5,2.5) +
	xlim(-2.5,2.5) +
	ggtitle("Accessibility vs Expression divergence - integrated categories")

	Z <- all_classes_integrated[all_classes_integrated$class == "inter",] %>%
	ggplot(aes(x=P_est.mean, y=P_est.mean_RNA, color=div_category)) +
	geom_point() +
	theme_main() +
	scale_color_discrete(guide = F) +
	xlab("Estimated accessibility divergence") +
	ylab("Estimated expression divergennce") +
	ylim(-2.5,2.5) +
	xlim(-2.5,2.5) +
	ggtitle("Accessibility vs Expression divergence - integrated categories")

	A <- all_classes_integrated[all_classes_integrated$class == "intra",] %>%
	ggplot(aes(x=P_est.mean, y=P_est.mean_RNA, color=div_category)) +
	geom_point() +
	theme_main() +
	scale_color_discrete(guide = F) +
	xlab("Estimated accessibility divergence") +
	ylab("Estimated expression divergennce") +
	ylim(-2.5,2.5) +
	xlim(-2.5,2.5) +
	ggtitle("Accessibility vs Expression divergence - integrated categories")

	B <- plot_grid(X, Y, Z, A)



		ggsave(G, file = "./AS_ATAC_RNA_2020_10_1/Figures_centered1000_runs/Acc_vs_Exp_categories_ZHR_Z30_noguide.pdf")

# % cis expression by group
M <- all_classes_integrated[all_classes_integrated$P_qvalue_RNA < 0.05 & all_classes_integrated$class == "start", ,] %>%
ggplot(aes(x=div_category, fill=div_category, y=perc_cis_RNA)) +
geom_boxplot(notch=TRUE) +
theme_main() +
ylab("Gene expression - percent cis") +
xlab("") +
scale_fill_discrete(guide=FALSE) +
scale_x_discrete(labels=c("Acc cons & Exp div","Acc div & Exp div\nopposite", "Acc div & Exp div\nsame")) +
ggtitle("Percent cis of expression divergence across integrated categories") +
theme(axis.text = element_text(size = 15),
axis.title = element_text(size = 15),
plot.title = element_text(size = 13, face = "bold"))
	ggsave(M, file = "./Figures/Exp_percCis_int_categories_ZHR_Z30.pdf")

# % cis accessibility by group
N <- all_classes_integrated[all_classes_integrated$P_qvalue < 0.05,] %>%
ggplot(aes(x=div_category, fill=div_category, y=perc_cis)) +
geom_boxplot(notch=TRUE) +
theme_main() +
ylab("Accessibility - percent cis") +
xlab("") +
scale_fill_discrete(guide=FALSE) +
scale_x_discrete(labels=c("Acc div & Exp cons","Acc div & Exp div\nopposite", "Acc div & Exp div\nsame")) +
ggtitle("Percent cis of accessibility across integrated categories") +
theme(axis.text = element_text(size = 15),
axis.title = element_text(size = 15),
plot.title = element_text(size = 13, face = "bold"))
	ggsave(N, file = "./Figures/Acc_percCis_int_categories_ZHR_Z30.pdf")

# effect size expression by group
L <- all_classes_integrated[all_classes_integrated$P_qvalue_RNA < 0.05,] %>%
ggplot(aes(x=div_category, fill=div_category, y=abs(P_est.mean_RNA))) +
geom_boxplot(notch=TRUE) +
ylim(0,1.5) +
ylab("Magnitude of estimated expression divergence") +
xlab("") +
scale_x_discrete(labels=c("Acc cons & Exp div","Acc div & Exp div\nopposite", "Acc div & Exp div\nsame")) +
scale_fill_discrete(guide=FALSE) +
theme_main() +
theme(axis.text = element_text(size = 15),
axis.title = element_text(size = 15),
plot.title = element_text(size = 13, face = "bold")) +
ggtitle("Magnitude of expression divergence across integrated categories")
	ggsave(L, file = "./Figures/Exp_effectsize_int_categories_ZHR_Z30.pdf")

# same as above but with JUST txStart to demonstrate that this isn't just due to incorrect pairing
K <- all_classes_integrated[all_classes_integrated$P_qvalue_RNA < 0.05 & all_classes_integrated$class == "start",] %>%
ggplot(aes(x=div_category, fill=div_category, y=abs(P_est.mean_RNA))) +
geom_boxplot(notch=TRUE) +
ylim(0,1.5) +
ylab("Magnitude of estimated expression divergence") +
xlab("") +
scale_x_discrete(labels=c("Acc cons & Exp div","Acc div & Exp div\nopposite", "Acc div & Exp div\nsame")) +
scale_fill_discrete(guide=FALSE) +
theme_main() +
theme(axis.text = element_text(size = 15),
axis.title = element_text(size = 15),
plot.title = element_text(size = 12, face = "bold")) +
ggtitle("Magnitude of expression divergence (only txStart) across integrated categories")
	ggsave(K, file = "./Figures/Exp_effectsize_txStartONLY_int_categories_ZHR_Z30.pdf")

# effect size accessibility by group
Y <- all_classes_integrated[all_classes_integrated$P_qvalue < 0.05,] %>%
ggplot(aes(x=div_category, fill=div_category, y=abs(P_est.mean))) +
geom_boxplot(notch=TRUE) +
ylim(0,1.5) +
ylab("Magnitude of estimated accessibility divergence") +
xlab("") +
scale_x_discrete(labels=c("Acc div & Exp cons","Acc div & Exp div\nopposite", "Acc div & Exp div\nsame")) +
scale_fill_discrete(guide=FALSE) +
theme_main() +
theme(axis.text = element_text(size = 15),
axis.title = element_text(size = 15),
plot.title = element_text(size = 13, face = "bold")) +
ggtitle("Magnitude of accessibility divergence across integrated categories")
	ggsave(Y, file = "./Figures/Acc_effectsize_int_categories_ZHR_Z30.pdf")

# effect size accessibility by group with ONLY txSTart
Z <- all_classes_integrated[all_classes_integrated$P_qvalue < 0.05 & all_classes_integrated$class == "start",] %>%
ggplot(aes(x=div_category, fill=div_category, y=abs(P_est.mean))) +
geom_boxplot(notch=TRUE) +
ylim(0,1.5) +
ylab("Magnitude of estimated accessibility divergence") +
xlab("") +
scale_x_discrete(labels=c("Acc div & Exp cons","Acc div & Exp div\nopposite", "Acc div & Exp div\nsame")) +
scale_fill_discrete(guide=FALSE) +
theme_main() +
theme(axis.text = element_text(size = 15),
axis.title = element_text(size = 15),
plot.title = element_text(size = 11, face = "bold")) +
ggtitle("Magnitude of accessibility divergence (only txStart) across integrated categories")
	ggsave(Z, file = "./Figures/Acc_effectsize_txStartONLY_int_categories_ZHR_Z30.pdf")

# compare distance to gene by category - limit to 30Kb
Q <- all_classes_integrated[(all_classes_integrated$class == "inter" | all_classes_integrated$class == "intra") & all_classes_integrated$div_category != "AccC_ExpC" & all_classes_integrated$div_category != "AccC_ExpD",] %>%
ggplot(aes(x=div_category, fill=div_category, y=distance_to_exon1/1000)) +
geom_boxplot(notch=TRUE) +
ylim(0,30) +
ylab("Distance to closest expressed gene") +
xlab("") +
#scale_x_discrete(labels=c("Acc div & Exp cons","Acc div & Exp div\nopposite", "Acc div & Exp div\nsame")) +
scale_fill_discrete(guide=FALSE) +
scale_x_discrete(labels=c("Acc div & Exp cons","Acc div & Exp div\nopposite", "Acc div & Exp div\nsame")) +
theme_main() +
theme(axis.text = element_text(size = 13),
axis.title = element_text(size = 15),
plot.title = element_text(size = 12, face = "bold")) +
ggtitle("Kilobases to closest gene (inter- & intragenic) across integrated categories")
	ggsave(Q, file = "./Figures/Dist_to_closest_gene_int_categories_ZHR_Z30.pdf")


class_props <- matrix(ncol = 4, nrow = 5) %>% as.data.frame()
colnames(class_props) <- c("inter", "intra", "start", "end")
class_props[1,1] <- nrow(all_classes_integrated[all_classes_integrated$class == "inter" & all_classes_integrated$div_category == "AccC_ExpC",])
class_props[2,1] <- nrow(all_classes_integrated[all_classes_integrated$class == "inter" & all_classes_integrated$div_category == "AccC_ExpD",])
class_props[3,1] <- nrow(all_classes_integrated[all_classes_integrated$class == "inter" & all_classes_integrated$div_category == "AccD_ExpC",])
class_props[4,1] <- nrow(all_classes_integrated[all_classes_integrated$class == "inter" & all_classes_integrated$div_category == "AccD_ExpD_O",])
class_props[5,1] <- nrow(all_classes_integrated[all_classes_integrated$class == "inter" & all_classes_integrated$div_category == "AccD_ExpD_S",])

class_props[1,2] <- nrow(all_classes_integrated[all_classes_integrated$class == "intra" & all_classes_integrated$div_category == "AccC_ExpC",])
class_props[2,2] <- nrow(all_classes_integrated[all_classes_integrated$class == "intra" & all_classes_integrated$div_category == "AccC_ExpD",])
class_props[3,2] <- nrow(all_classes_integrated[all_classes_integrated$class == "intra" & all_classes_integrated$div_category == "AccD_ExpC",])
class_props[4,2] <- nrow(all_classes_integrated[all_classes_integrated$class == "intra" & all_classes_integrated$div_category == "AccD_ExpD_O",])
class_props[5,2] <- nrow(all_classes_integrated[all_classes_integrated$class == "intra" & all_classes_integrated$div_category == "AccD_ExpD_S",])

class_props[1,3] <- nrow(all_classes_integrated[all_classes_integrated$class == "start" & all_classes_integrated$div_category == "AccC_ExpC",])
class_props[2,3] <- nrow(all_classes_integrated[all_classes_integrated$class == "start" & all_classes_integrated$div_category == "AccC_ExpD",])
class_props[3,3] <- nrow(all_classes_integrated[all_classes_integrated$class == "start" & all_classes_integrated$div_category == "AccD_ExpC",])
class_props[4,3] <- nrow(all_classes_integrated[all_classes_integrated$class == "start" & all_classes_integrated$div_category == "AccD_ExpD_O",])
class_props[5,3] <- nrow(all_classes_integrated[all_classes_integrated$class == "start" & all_classes_integrated$div_category == "AccD_ExpD_S",])

class_props[1,4] <- nrow(all_classes_integrated[all_classes_integrated$class == "end" & all_classes_integrated$div_category == "AccC_ExpC",])
class_props[2,4] <- nrow(all_classes_integrated[all_classes_integrated$class == "end" & all_classes_integrated$div_category == "AccC_ExpD",])
class_props[3,4] <- nrow(all_classes_integrated[all_classes_integrated$class == "end" & all_classes_integrated$div_category == "AccD_ExpC",])
class_props[4,4] <- nrow(all_classes_integrated[all_classes_integrated$class == "end" & all_classes_integrated$div_category == "AccD_ExpD_O",])
class_props[5,4] <- nrow(all_classes_integrated[all_classes_integrated$class == "end" & all_classes_integrated$div_category == "AccD_ExpD_S",])

class_props$inter <- class_props$inter/sum(class_props$inter)
class_props$intra <- class_props$intra/sum(class_props$intra)
class_props$start <- class_props$start/sum(class_props$start)
class_props$end <- class_props$end/sum(class_props$end)
class_props$div_category <- c("AccC_ExpC", "AccC_ExpD", "AccD_ExpC", "AccD_ExpD_O", "AccD_ExpD_S")

G <- melt(class_props) %>%
ggplot(aes(x=variable, y=value, fill=div_category)) +
geom_bar(stat="identity", colour="white") +
xlab("") +
ylab("Proportion") +
scale_x_discrete(labels=c("Intergenic","Intragenic", "txStart", "txEnd")) +
scale_fill_discrete(name = "Divergence category",
labels = c("Acc cons & Exp cons", "Acc cons & Exp div", "Acc div & Exp cons", "Acc div & Exp div opposite", "Acc div & Exp div same ")) +
theme_main() +
theme(axis.text = element_text(size = 12),
axis.title = element_text(size = 15),
plot.title = element_text(size = 13, face = "bold"),
legend.text= element_text(size = 10),
legend.title = element_text(size = 10)) +
ggtitle("Relative proportions of integrated category for each class")
	ggsave(G, file = "./Figures/Prop_div_categories_by_int_categories_ZHR_Z30.pdf")
