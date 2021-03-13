#!/bin/bash

# define files
strain1=ZHR
strain2=TSIM
hybrid=ZHR_TSIM
datatype=ATAC
ref_genome2="dm6"
final_peak_file1=/nfs/turbo/lsa-wittkopp/Lab/Henry/BEDfiles/dm6_txStart.bed
final_peak_file2=/nfs/turbo/lsa-wittkopp/Lab/Henry/BEDfiles/dm6_txEnd.bed
final_peak_file3=/nfs/turbo/lsa-wittkopp/Lab/Henry/BEDfiles/ZHR_Z30_TSIM_analysis_regions_intragenic_NO_TxStart_End_center1000.bed
final_peak_file4=/nfs/turbo/lsa-wittkopp/Lab/Henry/BEDfiles/ZHR_Z30_TSIM_analysis_regions_intergenic_NO_TxStart_End_center1000.bed

# move to working directory
cd /nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/MERGED_AND_FINAL_FILES/${datatype}

# Get read count for each sample/replicate for master peak set -- txStart
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/./multiBamCov -bams ../../ALL_${datatype}_ALIGNED/${strain1}_1_${datatype}/${strain1}_1_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_2_${datatype}/${strain1}_2_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_3_${datatype}/${strain1}_3_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_1_${datatype}/${strain2}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_2_${datatype}/${strain2}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_3_${datatype}/${strain2}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
-bed $final_peak_file1 \
> ${hybrid}_${datatype}_txStart500_counts.bed

# Get read count for each sample/replicate for master peak set -- txEND
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/./multiBamCov -bams ../../ALL_${datatype}_ALIGNED/${strain1}_1_${datatype}/${strain1}_1_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_2_${datatype}/${strain1}_2_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_3_${datatype}/${strain1}_3_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_1_${datatype}/${strain2}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_2_${datatype}/${strain2}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_3_${datatype}/${strain2}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
-bed $final_peak_file2 \
> ${hybrid}_${datatype}_txEnd500_counts.bed

# Get read count for each sample/replicate for master peak set -- intragenic
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/./multiBamCov -bams ../../ALL_${datatype}_ALIGNED/${strain1}_1_${datatype}/${strain1}_1_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_2_${datatype}/${strain1}_2_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_3_${datatype}/${strain1}_3_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_1_${datatype}/${strain2}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_2_${datatype}/${strain2}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_3_${datatype}/${strain2}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
-bed $final_peak_file3 \
> ${hybrid}_${datatype}_intragenic_peak_counts_center1000.bed

# Get read count for each sample/replicate for master peak set -- intergenic
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/./multiBamCov -bams ../../ALL_${datatype}_ALIGNED/${strain1}_1_${datatype}/${strain1}_1_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_2_${datatype}/${strain1}_2_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_3_${datatype}/${strain1}_3_${datatype}_tsim.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_1_${datatype}/${strain2}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_2_${datatype}/${strain2}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_3_${datatype}/${strain2}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
-bed $final_peak_file4 \
> ${hybrid}_${datatype}_intergenic_peak_counts_center1000.bed
