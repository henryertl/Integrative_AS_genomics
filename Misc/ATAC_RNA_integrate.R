# construct merged files from Bayes test outputs


# read in data
ATAC <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/Parental_test_output_ZHR_Z30_sub_ATAC_CPM_macs2_20min.txt", header = T, sep = " ")
RNA <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Bayes_outputs/Parental_test_output_RNA_full_CPM_20_1000.txt", header = T, sep = " ")
Exons <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/constitutive_regions_final_dm6.bed", header = F, sep = "\t")
