#!/bin/bash

#######################################################################################################################
#######################  GOAL : Align reads to concatenated genome, AS sort, and get BAM files ########################
#######################################################################################################################

# Run this script for each sample replicate in an allele specific NGS experiment
# Ensure directory structure is consistent with script
# Run as PART 1 before subsequent sequencing-specific scripts
#

######################## DEFINE INPUTS ##################################

strain1=
strain2=
strain1_chrom_names= # format e.g.: "chrom2L chrom 2R ..." -- ensure matches genome chromosome names
strain2_chrom_names=
#####
# this script is written to convert to a final genome in TWO parts (i.e. with an intermediate genome)
# adjust inputs and script as necessary if only need one (or more than one) conversion step(s)
ref_genome1=
ref_genome2=
strain1_to_refgenome1_chain=
strain2_to_refgenome1_chain=
refgenome1_to_refgenome2_chain=
refgenome2_fasta=
#####
indexed_cat_genome= #path/to/indexed_concatenated_genome
forward_reads= #/path/to/reads
reverse_reads= #/path/to/reads
SAMPLE_NAME= # e.g. ZHR_1
final_folder_output= # ensure this  directory exists with correct path structure

############## Change directory to sample specific folder #############
cd /nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/ALIGNED_ASsorted/$SAMPLE_NAME #ensure this directory exists

############# Align reads, clean up, and sort alleles ################
# Align reads
bowtie2 --very-sensitive -x $indexed_cat_genome \
        -1 $forward_reads \
        -2 $reverse_reads \
        -S ${SAMPLE_NAME}_${strain1}${strain2}cat.sam
samtools sort -O BAM -o ${SAMPLE_NAME}_${strain1}${strain2}cat.bam \
        ${SAMPLE_NAME}_${strain1}${strain2}cat.sam
samtools index ${SAMPLE_NAME}_${strain1}${strain2}cat.bam

## Clean up and grab Strain 1 chromosome reads
samtools view -b -h -f 3 -F 4 -F 256 -F 1024 -F 2048 -q 20 \
        ${SAMPLE_NAME}_${strain1}${strain2}cat.bam \
        $strain1_chrom_names \
        > ${SAMPLE_NAME}.${strain1}_sorted.bam

## Clean up and grab Strain 2 chromosome reads
samtools view -b -h -f 3 -F 4 -F 256 -F 1024 -F 2048 -q 20 \
        ${SAMPLE_NAME}_${strain1}${strain2}cat.bam \
        $strain2_chrom_names \
        > ${SAMPLE_NAME}.${strain2}_sorted.bam

# Sort by name
samtools sort -n -o ${SAMPLE_NAME}.${strain1}_sorted_name.bam \
        ${SAMPLE_NAME}.${strain1}_sorted.bam
samtools sort -n -o ${SAMPLE_NAME}.${strain2}_sorted_name.bam \
        ${SAMPLE_NAME}.${strain2}_sorted.bam

# Sort by coordinate
samtools sort -O BAM -o ${SAMPLE_NAME}.${strain1}_sorted_coord.bam \
        ${SAMPLE_NAME}.${strain1}_sorted.bam
mv ${SAMPLE_NAME}.${strain1}_sorted_coord.bam $final_folder_output

samtools sort -O BAM -o ${SAMPLE_NAME}.${strain2}_sorted_coord.bam \
        ${SAMPLE_NAME}.${strain2}_sorted.bam
mv ${SAMPLE_NAME}.${strain2}_sorted_coord.bam $final_folder_output

################  Convert BAM to BED with bedtools ##################
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/bamToBed -cigar -i ${SAMPLE_NAME}.${strain1}_sorted_name.bam \
> ${SAMPLE_NAME}.${strain1}_sorted_name.bed

/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/bamToBed -cigar -i ${SAMPLE_NAME}.${strain2}_sorted_name.bam \
> ${SAMPLE_NAME}.${strain2}_sorted_name.bed

################# Convert to reference genome #######################
# Strain 1 -> Reference genome 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${SAMPLE_NAME}.${strain1}_sorted_name.bed \
$strain1_to_refgenome1_chain \
${SAMPLE_NAME}.${strain1}_${ref_genome1}.bed \
${SAMPLE_NAME}.${strain1}_${ref_genome1}_unmapped

# Strain 2 -> Reference genome 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${SAMPLE_NAME}.${strain2}_sorted_name.bed \
$strain2_to_refgenome1_chain \
${SAMPLE_NAME}.${strain2}_${ref_genome1}.bed \
${SAMPLE_NAME}.${strain2}_${ref_genome1}_unmapped

# Reference genome 1 -> Reference genome 2 (Strain 1)
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${SAMPLE_NAME}.${strain1}_${ref_genome1}.bed \
$refgenome1_to_refgenome2_chain \
${SAMPLE_NAME}.${strain1}_${ref_genome2}.bed \
${SAMPLE_NAME}.${strain1}_${ref_genome2}_unmapped

# Reference genome 1 -> Reference genome 2 (Strain 2)
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${SAMPLE_NAME}.${strain2}_${ref_genome1}.bed \
$refgenome1_to_refgenome2_chain \
${SAMPLE_NAME}.${strain2}_${ref_genome2}.bed \
${SAMPLE_NAME}.${strain2}_${ref_genome2}_unmapped

#################### Get final BAM files #############################

# Convert BED to BAM
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/bedToBam -i ${SAMPLE_NAME}.${strain1}_${ref_genome2}.bed \
-g $refgenome2_fasta \
> ${SAMPLE_NAME}.${strain1}_${ref_genome2}.bam

/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/bedToBam -i ${SAMPLE_NAME}.${strain2}_${ref_genome2}.bed \
-g $refgenome2_fasta \
> ${SAMPLE_NAME}.${strain2}_${ref_genome2}.bam

# Sort and index final files
samtools sort -o ${SAMPLE_NAME}.${strain1}_${ref_genome2}_sorted_coord.bam \
        ${SAMPLE_NAME}.${strain1}_${ref_genome2}.bam
samtools index ${SAMPLE_NAME}.${strain1}_${ref_genome2}_sorted_coord.bam

samtools sort -o ${SAMPLE_NAME}.${strain2}_${ref_genome2}_sorted_coord.bam \
        ${SAMPLE_NAME}.${strain2}_${ref_genome2}.bam
samtools index ${SAMPLE_NAME}.${strain2}_${ref_genome2}_sorted_coord.bam
