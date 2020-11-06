#!/bin/bash
#PBS -N get_counts_genic_regions_allspp
#PBS -l nodes=2:ppn=2,mem=24gb,walltime=2:00:00
#PBS -S /bin/sh
#PBS -M hertl@umich.edu
#PBS -A lsa_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -j oe
#PBS -V


#### Define file inputs ####

# gene_beds_by_chrom #
chrX_genes_dm3=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/bed_files/chrX_genes_dm3_clean.bed
chr2L_genes_dm3=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/bed_files/chr2L_genes_dm3_clean.bed
chr2R_genes_dm3=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/bed_files/chr2R_genes_dm3_clean.bed
chr3R_genes_dm3=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/bed_files/chr3R_genes_dm3_clean.bed
chr3L_genes_dm3=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/bed_files/chr3L_genes_dm3_clean.bed
chr4_genes_dm3=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/bed_files/chr4_genes_dm3_clean.bed

# all genes all chroms
all_chroms_dm3=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/bed_files/all_genes_dm3_final_clean.bed

# bed genome files for each spp #
D_wil_dm3_bed=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_1_dm3_md_hqaa_sorted.bed
D_sim_dm3_bed=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_1_dm3_md_hqaa_sorted.bed
D_sech_dm3_bed=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_1_dm3_md_hqaa_sorted.bed
D_per_dm3_bed=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_1_dm3_md_hqaa_sorted.bed
D_vir_dm3_bed=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_1_dm3_md_hqaa_sorted.bed
D_yak_dm3_bed=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_1_dm3_md_hqaa_sorted.bed
D_pse_dm3_bed=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_1_dm3_md_hqaa_sorted.bed
D_ere_dm3_bed=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_1_dm3_md_hqaa_sorted.bed
D_moj_dm3_bed=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_1_dm3_md_hqaa_sorted.bed
D_ana_dm3_bed=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_1_dm3_md_hqaa_sorted.bed

# file outputs
D_wil_gene_counts_output=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_1_dm3_genes.bed
D_wil_gene_counts_output_clean=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_1_dm3_genes_final.bed
D_sim_gene_counts_output=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_1_dm3_genes.bed
D_sim_gene_counts_output_clean=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_1_dm3_genes_final.bed
D_sech_gene_counts_output=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_1_dm3_genes.bed
D_sech_gene_counts_output_clean=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_1_dm3_genes_final.bed
D_per_gene_counts_output=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_1_dm3_genes.bed
D_per_gene_counts_output_clean=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_1_dm3_genes_final.bed
D_vir_gene_counts_output=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_1_dm3_genes.bed
D_vir_gene_counts_output_clean=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_1_dm3_genes_final.bed
D_yak_gene_counts_output=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_1_dm3_genes.bed
D_yak_gene_counts_output_clean=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_1_dm3_genes_final.bed
D_pse_gene_counts_output=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_1_dm3_genes.bed
D_pse_gene_counts_output_clean=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_1_dm3_genes_final.bed
D_ere_gene_counts_output=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_1_dm3_genes.bed
D_ere_gene_counts_output_clean=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_1_dm3_genes_final.bed
D_moj_gene_counts_output=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_1_dm3_genes.bed
D_moj_gene_counts_output_clean=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_1_dm3_genes_final.bed
D_ana_gene_counts_output=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_1_dm3_genes.bed
D_ana_gene_counts_output_clean=/scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_1_dm3_genes_final.bed


if [ -e “$PBS_NODEFILE” ] ; then
   echo “Running on”
   uniq -c $PBS_NODEFILE
fi

if [ -d “$PBS_O_WORKDIR” ] ; then
   cd $PBS_O_WORKDIR
fi
echo “Running from $(pwd)”

# get counts from genic regions and clean up
##D_wil
bedmap --faster --echo --count $D_wil_dm3_bed $all_chroms_dm3 > $D_wil_gene_counts_output
sed 's/|/\t/g' $D_wil_gene_counts_output > $D_wil_gene_counts_output_clean
rm $D_wil_gene_counts_output

##D_per
bedmap --faster --echo --count $D_per_dm3_bed $all_chroms_dm3 > $D_per_gene_counts_output
sed 's/|/\t/g' $D_per_gene_counts_output > $D_per_gene_counts_output_clean
rm $D_per_gene_counts_output

##D_sim
bedmap --faster --echo --count $D_sim_dm3_bed $all_chroms_dm3 > $D_sim_gene_counts_output
sed 's/|/\t/g' $D_sim_gene_counts_output > $D_sim_gene_counts_output_clean
rm $D_sim_gene_counts_output

##D_vir
bedmap --faster --echo --count $D_sim_dm3_bed $all_chroms_dm3 > $D_vir_gene_counts_output
sed 's/|/\t/g' $D_sim_gene_counts_output > $D_vir_gene_counts_output_clean
rm $D_vir_gene_counts_output

##D_yak
bedmap --faster --echo --count $D_yak_dm3_bed $all_chroms_dm3 > $D_yak_gene_counts_output
sed 's/|/\t/g' $D_yak_gene_counts_output > $D_yak_gene_counts_output_clean
rm $D_yak_gene_counts_output

##D_pse
bedmap --faster --echo --count $D_yak_dm3_bed $all_chroms_dm3 > $D_pse_gene_counts_output
sed 's/|/\t/g' $D_yak_gene_counts_output > $D_pse_gene_counts_output_clean
rm $D_pse_gene_counts_output

##D_ere
bedmap --faster --echo --count $D_yak_dm3_bed $all_chroms_dm3 > $D_ere_gene_counts_output
sed 's/|/\t/g' $D_yak_gene_counts_output > $D_ere_gene_counts_output_clean
rm $D_ere_gene_counts_output

##D_moj
bedmap --faster --echo --count $D_moj_dm3_bed $all_chroms_dm3 > $D_moj_gene_counts_output
sed 's/|/\t/g' $D_moj_gene_counts_output > $D_moj_gene_counts_output_clean
rm $D_moj_gene_counts_output

##D_sech
bedmap --faster --echo --count $D_sech_dm3_bed $all_chroms_dm3 > $D_sech_gene_counts_output
sed 's/|/\t/g' $D_sech_gene_counts_output > $D_sech_gene_counts_output_clean
rm $D_sech_gene_counts_output

##D_ana
bedmap --faster --echo --count $D_ana_dm3_bed $all_chroms_dm3 > $D_ana_gene_counts_output
sed 's/|/\t/g' $D_ana_gene_counts_output > $D_ana_gene_counts_output_clean
rm $D_ana_gene_counts_output
