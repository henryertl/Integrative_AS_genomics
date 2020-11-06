#!/bin/bash

# numerous diff SNP files and programs need to be used to generate the SNP/variation files needed for ASE analysis
## more complicated than I thought it would be... so going to keep a record here of each step

# 1) align gDNA of sp 1 to sp2 genome with bowtie default parameters


# 2) call SNPs with gatk haplotypeCaller


# 3) split up vcf file from output file into diff chroms
for i in { zhr_bow_chr2L zhr_bow_chr2R zhr_bow_chr3L zhr_bow_chr3R zhr_bow_chr4 zhr_bow_chrX };
do vcftools  --vcf  tsim_zhr_variants_longer.vcf  --chr $i  --recode --recode-INFO-all --out  VCF_$i;
done

# 4) gzip all output files
gzip -r *.vcf /path/to/directory/with_vcf_files

#5) create text based SNP file with WASP script
## this file is good for downstream WASP applications
./extract_vcf_snps.sh /scratch/lsa_root/lsa/hertl/genomes/SNP_files/test_run /scratch/lsa_root/lsa/hertl/genomes/SNP_files/test_run ##ensure that names are edited to chromosome names from masked genome file for downstream WASP applications 

## to generate alteration for ASEr applications, edit file(s) with following code
gunzip *.snps.txt.gz
awk '{print $1"\t"$1+1"\t"$2"|"$3}' chr2L.snps.txt > chr2L.snps.bed
awk '{print $1"\t"$1+1"\t"$2"|"$3}' chr2R.snps.txt > chr2R.snps.bed
awk '{print $1"\t"$1+1"\t"$2"|"$3}' chr3L.snps.txt > chr3L.snps.bed
awk '{print $1"\t"$1+1"\t"$2"|"$3}' chr3R.snps.txt > chr3R.snps.bed
awk '{print $1"\t"$1+1"\t"$2"|"$3}' chr4.snps.txt > chr4.snps.bed
awk '{print $1"\t"$1+1"\t"$2"|"$3}' chrX.snps.txt > chrX.snps.bed
## to generate same files as above but with chrom names included, and then also one combined file
awk '{print "zhr_bow_chr2L""\t"$1"\t"$1+1"\t"$2"|"$3}' chr2L.snps.txt > chr2L_names.snps.bed
awk '{print "zhr_bow_chr2R""\t"$1"\t"$1+1"\t"$2"|"$3}' chr2R.snps.txt > chr2R_names.snps.bed
awk '{print "zhr_bow_chr3L""\t"$1"\t"$1+1"\t"$2"|"$3}' chr3L.snps.txt > chr3L_names.snps.bed
awk '{print "zhr_bow_chr3R""\t"$1"\t"$1+1"\t"$2"|"$3}' chr3R.snps.txt > chr3R_names.snps.bed
awk '{print "zhr_bow_chr4""\t"$1"\t"$1+1"\t"$2"|"$3}' chr4.snps.txt > chr4_names.snps.bed
awk '{print "zhr_bow_chrX""\t"$1"\t"$1+1"\t"$2"|"$3}' chrX.snps.txt > chrX_names.snps.bed
cat chr2L_names.snps.bed chr2R_names.snps.bed chr3L_names.snps.bed chr3R_names.snps.bed chr4_names.snps.bed chrX_names.snps.bed > chrALL_names.snps.bed

# generate SNP_split SNP file
awk '{print "1""\t"$1"\t"$2"\t""1""\t"$4}' chrALL_names.snps.bed | tr '|' '/' > chrALL_names.snpsplit.txt

# SNP_split command
./SNPsplit --snp_file /scratch/lsa_root/lsa/hertl/genomes/SNP_files/SNPsplit_files/chrALL_names.snpsplit.txt /scratch/lsa_root/lsa/hertl/executables/WASP/map1/zhr_z30_temponmode.bam --paired --no_sort


/scratch/lsa_root/lsa/hertl/other_datasets/tempo_mode_gDNA/TSIM_ZHR_hyb_peRNAseq_F.fq
