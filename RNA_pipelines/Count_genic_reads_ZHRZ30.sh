#!/bin/bash

#########################################################################################
#####################  GOAL : Count genic reads and normalize  ##########################
#########################################################################################


################################ Define files ####################################

strain1=ZHR
strain2=Z30
hybrid=ZHR_Z30
exon_coordiantes_BED_file=/nfs/turbo/lsa-wittkopp/Lab/Henry/BEDfiles/exons_dm6.bed
ref_genome2="dm6"
datatype="RNA"

# Move to working directory
cd /nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/MERGED_AND_FINAL_FILES/${datatype}

############### Get table with genic read counts for each sample #################
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/./multiBamCov -bams ../../ALL_${datatype}_ALIGNED/${strain1}_1_${datatype}/${strain1}_1_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_2_${datatype}/${strain1}_2_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_3_${datatype}/${strain1}_3_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_4_${datatype}/${strain1}_4_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_5_${datatype}/${strain1}_5_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_6_${datatype}/${strain1}_6_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_1_${datatype}/${strain2}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_2_${datatype}/${strain2}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_3_${datatype}/${strain2}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_4_${datatype}/${strain2}_4_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_5_${datatype}/${strain2}_5_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_6_${datatype}/${strain2}_6_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_4_${datatype}/${hybrid}_4_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_5_${datatype}/${hybrid}_5_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_6_${datatype}/${hybrid}_6_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_4_${datatype}/${hybrid}_4_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_5_${datatype}/${hybrid}_5_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_6_${datatype}/${hybrid}_6_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
-bed $exon_coordiantes_BED_file > ${hybrid}_exon_counts_${ref_genome2}.txt

# R script to add header and sum up exons
Rscript /path/to/Exon_sum.R
