# Generate count file for constitutive exons
strain1=ZHR
strain2=Z30
hybrid=ZHR_Z30
exon_coordiantes_BED_file=/nfs/turbo/lsa-wittkopp/Lab/Henry/BEDfiles/constitutive_regions_final_dm6.bed
ref_genome2="dm6"

cd /nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project

############### Get table with genic read counts for each sample #################
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/./multiBamCov -bams ALL_RNA_ALIGNED/${strain1}_1_RNA/${strain1}_1_RNA.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${strain1}_3_RNA/${strain1}_3_RNA.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${strain1}_5_RNA/${strain1}_5_RNA.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${strain2}_1_RNA/${strain2}_1_RNA.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${strain2}_3_RNA/${strain2}_3_RNA.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${strain2}_5_RNA/${strain2}_5_RNA.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${hybrid}_1_RNA/${hybrid}_1_RNA.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${hybrid}_3_RNA/${hybrid}_3_RNA.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${hybrid}_5_RNA/${hybrid}_5_RNA.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${hybrid}_1_RNA/${hybrid}_1_RNA.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${hybrid}_3_RNA/${hybrid}_3_RNA.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_RNA_ALIGNED/${hybrid}_5_RNA/${hybrid}_5_RNA.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${strain1}_1_ATAC/${strain1}_1_ATAC.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${strain1}_2_ATAC/${strain1}_2_ATAC.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${strain1}_3_ATAC/${strain1}_3_ATAC.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${strain2}_1_ATAC/${strain2}_1_ATAC.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${strain2}_2_ATAC/${strain2}_2_ATAC.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${strain2}_3_ATAC/${strain2}_3_ATAC.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${hybrid}_1_ATAC/${hybrid}_1_ATAC.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${hybrid}_2_ATAC/${hybrid}_2_ATAC.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${hybrid}_3_ATAC/${hybrid}_3_ATAC.${strain1}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${hybrid}_1_ATAC/${hybrid}_1_ATAC.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${hybrid}_2_ATAC/${hybrid}_2_ATAC.${strain2}_${ref_genome2}_sorted_coord.bam \
ALL_ATAC_ALIGNED/${hybrid}_3_ATAC/${hybrid}_3_ATAC.${strain2}_${ref_genome2}_sorted_coord.bam \
-bed $exon_coordiantes_BED_file > MERGED_AND_FINAL_FILES/ATAC_RNA_integrated_files/${hybrid}_exon_counts_ATAC_RNA_${ref_genome2}.txt

# Run Rscript to get remove < 20 (and >1000 for RNA) and compute CPM across all
Rscript /nfs/turbo/lsa-wittkopp/Lab/Henry/scripts/ATAC_RNA_CPM_transform.R
