##Libraries
library(data.table)
library(plyr)
library(dplyr)
library(INLA)
library(parallel)
library(qvalue)
library(magrittr)


##############################
##Get command line arguments##
##############################
Arguments = (commandArgs(TRUE))

#########################
##Define test functions##
#########################
##################
##Parental model##
##################
Parental_model <- function(locus_ID, Parental_data){

	##ADJUST WHEN WE HAVE FULL DATA
	chrom <- Parental_data[Parental_data$Paste_locus == locus_ID,]$chrom
	pos_start <- Parental_data[Parental_data$Paste_locus == locus_ID,]$start
	pos_end <- Parental_data[Parental_data$Paste_locus == locus_ID,]$end
	Paste_locus <- Parental_data$Paste_locus[Parental_data$Paste_locus == locus_ID]

	##Reformat data for modelling
	reformed_matrix <- matrix(ncol = 2, nrow = 6) %>% as.data.frame()
	reformed_matrix[1,] <- c(Parental_data[Parental_data$Paste_locus == locus_ID, c(5, 11)])
	reformed_matrix[2,] <- c(Parental_data[Parental_data$Paste_locus == locus_ID, c(6, 12)])
	reformed_matrix[3,] <- c(Parental_data[Parental_data$Paste_locus == locus_ID, c(7, 13)])
	reformed_matrix[4,] <- c(Parental_data[Parental_data$Paste_locus == locus_ID, c(8, 14)])
	reformed_matrix[5,] <- c(Parental_data[Parental_data$Paste_locus == locus_ID, c(9, 15)])
	reformed_matrix[6,] <- c(Parental_data[Parental_data$Paste_locus == locus_ID, c(10, 16)])

	names(reformed_matrix)[1] <- "p1_reads"
	names(reformed_matrix)[2] <- "p2_reads"
	reformed_matrix$Total_reads <- reformed_matrix$p1_reads + reformed_matrix$p2_reads

	P_mod <- inla(p1_reads ~ 1, data = reformed_matrix , family = "binomial", Ntrials = Total_reads)

	coef <- P_mod$summary.fixed
	fixed_effect_posterior <- P_mod$marginals.fixed[[1]]
	lower_p <- inla.pmarginal(0, fixed_effect_posterior)
	upper_p <- 1 - inla.pmarginal(0, fixed_effect_posterior)
	post_pred_p <- 2 * (min(lower_p, upper_p))


	P_mod_output <- data.table(chrom = chrom, pos_start = pos_start, pos_end = pos_end, Paste_locus = Paste_locus, P_est = coef[1], P_p_value = post_pred_p)

	return(P_mod_output)
}

################
##Hybrid model##
################
Hybrid_model <- function(locus_ID, Hybrid_data){

	##ADJUST WHEN WE HAVE FULL DATA
	chrom <- Hybrid_data[Hybrid_data$Paste_locus == locus_ID,]$chrom
	pos_start <- Hybrid_data[Hybrid_data$Paste_locus == locus_ID,]$start
	pos_end <- Hybrid_data[Hybrid_data$Paste_locus == locus_ID,]$end
	Paste_locus <- Hybrid_data$Paste_locus[Hybrid_data$Paste_locus == locus_ID]

	##ADJUST WHEN WE HAVE FULL DATA
	##Reformat data for modelling
	reformed_matrix <- matrix(ncol = 2, nrow = 6) %>% as.data.frame()
	reformed_matrix[1,] <- c(Hybrid_data[Hybrid_data$Paste_locus == locus_ID, c(5, 11)])
	reformed_matrix[2,] <- c(Hybrid_data[Hybrid_data$Paste_locus == locus_ID, c(6, 12)])
	reformed_matrix[3,] <- c(Hybrid_data[Hybrid_data$Paste_locus == locus_ID, c(7, 13)])
	reformed_matrix[4,] <- c(Hybrid_data[Hybrid_data$Paste_locus == locus_ID, c(8, 14)])
	reformed_matrix[5,] <- c(Hybrid_data[Hybrid_data$Paste_locus == locus_ID, c(9, 15)])
	reformed_matrix[6,] <- c(Hybrid_data[Hybrid_data$Paste_locus == locus_ID, c(10, 16)])

	names(reformed_matrix)[1] <- "p1_reads"
	names(reformed_matrix)[2] <- "p2_reads"
	reformed_matrix$Total_reads <-  reformed_matrix$p1_reads + reformed_matrix$p2_reads

	H_mod <- inla(p1_reads ~ 1, data = reformed_matrix , family = "binomial", Ntrials = Total_reads)

	coef <- H_mod$summary.fixed
	fixed_effect_posterior <- H_mod$marginals.fixed[[1]]
	lower_p <- inla.pmarginal(0, fixed_effect_posterior)
	upper_p <- 1 - inla.pmarginal(0, fixed_effect_posterior)
	post_pred_p <- 2 * (min(lower_p, upper_p))


	H_mod_output <- data.table(chrom = chrom, pos_start = pos_start, pos_end = pos_end, Paste_locus = Paste_locus, H_est = coef[1], H_p_value = post_pred_p)

	return(H_mod_output)
}


###########################
##Parental - Hybrid model##
###########################
Parental_hybrid_model <- function(locus_ID, Parental_hybrid_data){

	##ADJUST WHEN WE HAVE FULL DATA
	chrom <- Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID,]$chrom
	pos_start <- Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID,]$start
	pos_end <- Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID,]$end
	Paste_locus <- Parental_hybrid_data$Paste_locus[Parental_hybrid_data$Paste_locus == locus_ID]

	##ADJUST WHEN WE HAVE FULL DATA
	##Reformat data for modelling
	reformed_matrix <- matrix(ncol = 2, nrow = 12) %>% as.data.frame()
	reformed_matrix[1,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(5, 11)]) #p1
	reformed_matrix[2,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(6, 12)]) #p2
	reformed_matrix[3,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(7, 13)])	#h1
	reformed_matrix[4,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(8, 14)]) #h2
	reformed_matrix[5,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(9, 15)]) #p1
	reformed_matrix[6,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(10, 16)]) #p2
	reformed_matrix[7,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(17, 23)])	#h1
	reformed_matrix[8,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(18, 24)]) #h2
	reformed_matrix[9,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(19, 25)]) #p1
	reformed_matrix[10,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(20, 26)]) #p2
	reformed_matrix[11,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(21, 27)])	#h1
	reformed_matrix[12,] <- c(Parental_hybrid_data[Parental_hybrid_data$Paste_locus == locus_ID, c(22, 28)]) #h2


	names(reformed_matrix)[1] <- "p1_reads"
	names(reformed_matrix)[2] <- "p2_reads"
	reformed_matrix$Total_reads <-  reformed_matrix$p1_reads + reformed_matrix$p2_reads
	reformed_matrix$Env <-  c("Parental", "Parental", "Parental", "Parental", "Parental", "Parental", "Hybrid", "Hybrid", "Hybrid", "Hybrid", "Hybrid", "Hybrid")

	H_P_mod <- inla(p1_reads ~ Env, data = reformed_matrix , family = "binomial", Ntrials = Total_reads)

	coef <- H_P_mod$summary.fixed
	fixed_effect_posterior <- H_P_mod$marginals.fixed[[2]] ##Check this... and run some tests
	lower_p <- inla.pmarginal(0, fixed_effect_posterior)
	upper_p <- 1 - inla.pmarginal(0, fixed_effect_posterior)
	post_pred_p <- 2 * (min(lower_p, upper_p))

	H_P_mod_output <- data.table(chrom = chrom, pos_start = pos_start, pos_end = pos_end, Paste_locus = Paste_locus, H_P_est = coef[2,1], H_P_p_value = post_pred_p)

	return(H_P_mod_output)
}

######################
##Read primary data ##
######################
full_dataset <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/ZHR_TSIM_genic_counts_CPM_final_dm6_20min.txt", header = T)
full_dataset$chrom <- "0"
full_dataset$start <- "0"
full_dataset$end <- "0"

full_dataset$Paste_locus <- paste(full_dataset$chrom, full_dataset$start, full_dataset$end, full_dataset$gene, sep = "_")

## Get sums
All_reads_p1_r1 <- sum(full_dataset$P1_1)
All_reads_p1_r2 <- sum(full_dataset$P1_2)
All_reads_p1_r3 <- sum(full_dataset$P1_3)
All_reads_p1_r4 <- sum(full_dataset$P1_4)
All_reads_p1_r5 <- sum(full_dataset$P1_5)
All_reads_p1_r6 <- sum(full_dataset$P1_6)

All_reads_p2_r1 <- sum(full_dataset$P2_1)
All_reads_p2_r2 <- sum(full_dataset$P2_2)
All_reads_p2_r3 <- sum(full_dataset$P2_3)
All_reads_p2_r4 <- sum(full_dataset$P2_4)
All_reads_p2_r5 <- sum(full_dataset$P2_5)
All_reads_p2_r6 <- sum(full_dataset$P2_6)

All_reads_p1_hyb_r1 <- sum(full_dataset$HYB_1_P1)
All_reads_p1_hyb_r2 <- sum(full_dataset$HYB_2_P1)
All_reads_p1_hyb_r3 <- sum(full_dataset$HYB_3_P1)
All_reads_p1_hyb_r4 <- sum(full_dataset$HYB_4_P1)
All_reads_p1_hyb_r5 <- sum(full_dataset$HYB_5_P1)
All_reads_p1_hyb_r6 <- sum(full_dataset$HYB_6_P1)
All_reads_p2_hyb_r1 <- sum(full_dataset$HYB_1_P2)
All_reads_p2_hyb_r2 <- sum(full_dataset$HYB_2_P2)
All_reads_p2_hyb_r3 <- sum(full_dataset$HYB_3_P2)
All_reads_p2_hyb_r4 <- sum(full_dataset$HYB_4_P2)
All_reads_p2_hyb_r5 <- sum(full_dataset$HYB_4_P2)
All_reads_p2_hyb_r6 <- sum(full_dataset$HYB_4_P2)

##Get collumns for Parent, hybrid and hybrid parent data

Parental_data <- full_dataset[, c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P1_4", "P1_5", "P1_6", "P2_1", "P2_2", "P2_3", "P2_4", "P2_5", "P2_6")]

Hybrid_data <- full_dataset[, c("chrom", "start", "end", "Paste_locus", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_4_P1", "HYB_5_P1", "HYB_6_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "HYB_4_P2", "HYB_5_P2", "HYB_6_P2")]

Parental_hybrid_data <- full_dataset[, c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P1_4", "P1_5", "P1_6", "P2_1", "P2_2", "P2_3", "P2_4", "P2_5", "P2_6",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_4_P1", "HYB_5_P1", "HYB_6_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "HYB_4_P2", "HYB_5_P2", "HYB_6_P2")]


#############################################
##run approrpriate test based on argument 1##
#############################################
if (Arguments[1] == "Parents") {

  Parental_results <- do.call(rbind, mclapply(Parental_data$Paste_locus, function(x) Parental_model(x, Parental_data), mc.cores = 4))

  write.table(Parental_results, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Bayes_outputs/Parental_test_output_ZHR_TSIM_RNA_full_CPM_20.txt", row.names = F, quote = F)

} else if (Arguments[1] == "Hybrids"){

  Hybrid_results <- do.call(rbind, mclapply(Hybrid_data$Paste_locus, function(x) Hybrid_model(x, Hybrid_data), mc.cores = 4))

  write.table(Hybrid_results, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Bayes_outputs/Hybrid_test_output_ZHR_TSIM_RNA_full_CPM_20.txt", row.names = F, quote = F)

} else if(Arguments[1] == "Parent-Hybrid") {

  Parental_hybrid_results <- do.call(rbind, mclapply(Parental_hybrid_data$Paste_locus, function(x) Parental_hybrid_model(x, Parental_hybrid_data), mc.cores = 4))

  write.table(Parental_hybrid_results, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Bayes_outputs/Parental_Hybrid_test_output_ZHR_TSIM_RNA_full_CPM_20.txt", row.names = F, quote = F)

}
