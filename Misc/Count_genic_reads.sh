#!/bin/bash

#########################################################################################
#####################  GOAL : Count genic reads and normalize  ##########################
#########################################################################################


################################ Define files ####################################

strain1=
strain2=
hybrid=
gene_coordiantes=
genome=

# Move to working directory
cd /scratch/lsa_root/lsa/hertl/AS_RNA_seq/ZHR_Z30_gene_counts

############### Get table with genic read counts for each sample #################

bedtools multicov -bams ${strain1}_1_.${strain1}_${genome}_sorted_coord.bam \
${strain1}_2_.${strain1}_${genome}_sorted_coord.bam \
${strain1}_3_.${strain1}_${genome}_sorted_coord.bam \
${strain1}_4_.${strain1}_${genome}_sorted_coord.bam \
${strain1}_5_.${strain1}_${genome}_sorted_coord.bam \
${strain1}_6_.${strain1}_${genome}_sorted_coord.bam \
${strain2}_1_.${strain1}_${genome}_sorted_coord.bam \
${strain2}_2_.${strain1}_${genome}_sorted_coord.bam \
${strain2}_3_.${strain1}_${genome}_sorted_coord.bam \
${strain2}_4_.${strain1}_${genome}_sorted_coord.bam \
${strain2}_5_.${strain1}_${genome}_sorted_coord.bam \
${strain2}_6_.${strain1}_${genome}_sorted_coord.bam \
${hybrid}_1_.${strain1}_${genome}_sorted_coord.bam \
${hybrid}_2_.${strain1}_${genome}_sorted_coord.bam \
${hybrid}_3_.${strain1}_${genome}_sorted_coord.bam \
${hybrid}_4_.${strain1}_${genome}_sorted_coord.bam \
${hybrid}_5_.${strain1}_${genome}_sorted_coord.bam \
${hybrid}_6_.${strain1}_${genome}_sorted_coord.bam \
${hybrid}_1_.${strain2}_${genome}_sorted_coord.bam \
${hybrid}_2_.${strain2}_${genome}_sorted_coord.bam \
${hybrid}_3_.${strain2}_${genome}_sorted_coord.bam \
${hybrid}_4_.${strain2}_${genome}_sorted_coord.bam \
${hybrid}_5_.${strain2}_${genome}_sorted_coord.bam \
${hybrid}_6_.${strain2}_${genome}_sorted_coord.bam \
-bed $gene_coordiantes > ${hybrid}_genic_counts_${genome}.bed

####################### Normalize via CPM (???) ############################
