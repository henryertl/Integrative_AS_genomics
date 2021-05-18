setwd("/Users/wittkopp_member/Code/")

# read in SNP files
Z30_NA_bed <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/pop_genome_analyses/Z30_gDNA_aligned_to_I01_ALL.bed", header = F)
ZHR_NA_bed <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/pop_genome_analyses/ZHR_gDNA_aligned_to_I01_ALL.bed", header = F)

Z30_ZIM_bed <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/pop_genome_analyses/Z30_gDNA_aligned_to_ZS10_ALL.bed", header = F)
ZHR_ZIM_bed <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/pop_genome_analyses/ZHR_gDNA_aligned_to_ZS10_ALL.bed", header = F)

snps <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/BED_files_for_analyses/zhr_z30_SNPs_dm6_genome.bed", header = F)
