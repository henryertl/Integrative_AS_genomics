#!/bin/bash

#######################################################################################################################
#######################  GOAL : Align reads to concatenated genome, AS sort, and get BAM files ########################
#######################################################################################################################

# Run this script for each sample replicate in an allele specific NGS experiment
# Ensure directory structure is consistent with script
# Run as PART 1 before subsequent sequencing-specific scripts
#

######################## DEFINE INPUTS ##################################

strain1="ZHR"
strain2="SIM"
strain1_chrom_names="zhr_bow_chrX zhr_bow_chr2L zhr_bow_chr2R zhr_bow_chr3L zhr_bow_chr3R zhr_bow_chr4"
strain2_chrom_names="sim_bow_chrX sim_bow_chr2L sim_bow_chr2R sim_bow_chr3L sim_bow_chr3R sim_bow_chr4"
#####
# this script is written to convert to a final genome in TWO parts (i.e. with an intermediate genome)
# adjust inputs and script as necessary if only need one (or more than one) conversion step(s)
ref_genome1="dm3"
ref_genome1B="sim1"
ref_genome2="dm6"
strain1_to_refgenome1_chain=/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/chainfiles/zhr_bow_to_dm3.chain
strain2_to_refgenome1B_chain=/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/chainfiles/sim_bow_to_drosim1.chain
refgenome1B_to_refgenome1_chain=/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/chainfiles/droSim1ToDm3.over.chain
refgenome1_to_refgenome2_chain=/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/chainfiles/dm3ToDm6.over.chain
refgenome2_genome_annotation=/nfs/turbo/lsa-wittkopp/Lab/Henry/genome_fasta_and_index/dm6.genome
#####
indexed_cat_genome=/nfs/turbo/lsa-wittkopp/Lab/Henry/genome_fasta_and_index/zhr_tsim_cat_bowtie_index #path/to/indexed_concatenated_genome
SAMPLE_NAME="ZHR_TSIM_5_RNA" # e.g. ZHR_1
forward_reads=/scratch/lsa_root/lsa/hertl/1777-HE/TRIMMED_READS/*S358_R1_001_val_1.fq.gz #/path/to/reads
reverse_reads=/scratch/lsa_root/lsa/hertl/1777-HE/TRIMMED_READS/*S358_R2_001_val_2.fq.gz #/path/to/reads

#final_folder_output= # ensure this  directory exists with correct path structure

############## Change directory to sample specific folder #############
cd /nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/ALL_RNA_ALIGNED/$SAMPLE_NAME #ensure this directory exists

############# Align reads, clean up, and sort alleles ################
# Align reads
bowtie2 -p 8 --very-sensitive-local -x $indexed_cat_genome \
        -1 $forward_reads \
        -2 $reverse_reads \
        -S ${SAMPLE_NAME}_${strain1}${strain2}cat.sam
samtools sort -O BAM -o ${SAMPLE_NAME}_${strain1}${strain2}cat.bam \
        ${SAMPLE_NAME}_${strain1}${strain2}cat.sam
samtools index ${SAMPLE_NAME}_${strain1}${strain2}cat.bam

## Clean up and grab Strain 1 chromosome reads
samtools view -b -h -f 3 -F 4 -F 256 -F 1024 -F 2048 -q 10 \
        ${SAMPLE_NAME}_${strain1}${strain2}cat.bam \
        $strain1_chrom_names \
        > ${SAMPLE_NAME}.${strain1}_sorted.bam

## Clean up and grab Strain 2 chromosome reads
samtools view -b -h -f 3 -F 4 -F 256 -F 1024 -F 2048 -q 10 \
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

samtools sort -O BAM -o ${SAMPLE_NAME}.${strain2}_sorted_coord.bam \
        ${SAMPLE_NAME}.${strain2}_sorted.bam

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

# Strain 2 -> Reference genome 1B
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${SAMPLE_NAME}.${strain2}_sorted_name.bed \
$strain2_to_refgenome1B_chain \
${SAMPLE_NAME}.${strain2}_${ref_genome1B}.bed \
${SAMPLE_NAME}.${strain2}_${ref_genome1B}_unmapped

## Reference genome 1B -> Reference genome 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${SAMPLE_NAME}.${strain2}_${ref_genome1B}.bed \
$refgenome1B_to_refgenome1_chain \
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
-g $refgenome2_genome_annotation \
> ${SAMPLE_NAME}.${strain1}_${ref_genome2}.bam

/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/bedToBam -i ${SAMPLE_NAME}.${strain2}_${ref_genome2}.bed \
-g $refgenome2_genome_annotation \
> ${SAMPLE_NAME}.${strain2}_${ref_genome2}.bam

# Sort and index final files
samtools sort -o ${SAMPLE_NAME}.${strain1}_${ref_genome2}_sorted_coord.bam \
        ${SAMPLE_NAME}.${strain1}_${ref_genome2}.bam
samtools index ${SAMPLE_NAME}.${strain1}_${ref_genome2}_sorted_coord.bam

samtools sort -o ${SAMPLE_NAME}.${strain2}_${ref_genome2}_sorted_coord.bam \
        ${SAMPLE_NAME}.${strain2}_${ref_genome2}.bam
samtools index ${SAMPLE_NAME}.${strain2}_${ref_genome2}_sorted_coord.bam
