#!/bin/bash

##############################################################################################
######### Goal: Construct final peak set and grab read counts for each replication  ##########
##############################################################################################

#################### Define files #####################
strain1=
strain2=
hybrid=
datatype=
ref_genome1=
ref_genome2=
strain1_to_refgenome1_chain=
strain2_to_refgenome1_chain=
refgenome1_to_refgenome2_chain=
genome2_genome_fasta=

# Move to working directory
cd /nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/MERGED_AND_FINAL_FILES/ATACseq

################# Merge replicates for each genotype  #########################
## Strain 1
samtools merge -f ${strain1}_merged_${strain1}allele.bam \
../${strain1}_1/${strain1}_1_${strain1}_sorted_coord.bam \
../${strain1}_2/${strain1}_2_${strain1}_sorted_coord.bam \
../${strain1}_3/${strain1}_3_${strain1}_sorted_coord.bam

samtools sort -o ${strain1}_merged_${strain1}allele.bam \
${strain1}_merged_${strain1}allele.bam
samtools index ${strain1}_merged_${strain1}allele.bam

## Strain 2
samtools merge -f ${strain2}_merged_${strain2}allele.bam \
../${strain2}_1/${strain2}_1_${strain2}_sorted_coord.bam \
../${strain2}_2/${strain2}_2_${strain2}_sorted_coord.bam \
../${strain2}_3/${strain2}_3_${strain2}_sorted_coord.bam

samtools sort -o ${strain2}_merged_${strain2}allele.bam \
${strain2}_merged_${strain2}allele.bam
samtools index ${strain2}_merged_${strain2}allele.bam

## Hybrid allele - Strain 1
samtools merge -f ${hybrid}_merged_${strain1}allele.bam \
../${hybrid}_1/${hybrid}_1_${strain1}_sorted_coord.bam \
../${hybrid}_2/${hybrid}_2_${strain1}_sorted_coord.bam \
../${hybrid}_3/${hybrid}_3_${strain1}_sorted_coord.bam

samtools sort -o ${hybrid}_merged_${strain1}allele.bam \
${hybrid}_merged_${strain1}allele.bam
samtools index ${hybrid}_merged_${strain1}allele.bam

## Hybrid allele - Strain 2
samtools merge -f ${hybrid}_merged_${strain2}allele.bam \
../${hybrid}_1/${hybrid}_1_${strain2}_sorted_coord.bam \
../${hybrid}_2/${hybrid}_2_${strain2}_sorted_coord.bam \
../${hybrid}_3${hybrid}_3_${strain2}_sorted_coord.bam

samtools sort -o ${hybrid}_merged_${strain2}allele.bam \
${hybrid}_merged_${strain2}allele.bam
samtools index ${hybrid}_merged_${strain2}allele.bam

############### Call peaks with HMMRATAC ##########################
## Generate genome files needed for HMMRATAC
samtools view -H ${strain1}_merged_${strain1}allele.bam \
| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' \
> ${strain1}.genome

samtools view -H ${strain2}_merged_${strain2}allele.bam \
| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' \
> ${strain2}.genome

## Run HMMRATAC
### Strain 1
java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/HMMRATAC_V1.2.5_exe.jar \
-b ${strain1}_merged_${strain1}allele.bam \
-i ${strain1}_merged_${strain1}allele.bam.bai \
-g ${strain1}.genome \
-o ${strain1}_merged_${strain1}allele.HMMRATAC \
--bedgraph TRUE

### Strain 2
java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/HMMRATAC_V1.2.5_exe.jar \
-b ${strain2}_merged_${strain2}allele.bam \
-i ${strain2}_merged_${strain2}allele.bam.bai \
-g ${strain2}.genome \
-o ${strain2}_merged_${strain2}allele.HMMRATAC \
--bedgraph TRUE

### Hybrid allele - Strain 1
java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/HMMRATAC_V1.2.5_exe.jar \
-b ${hybrid}_merged_${strain1}allele.bam \
-i ${hybrid}_merged_${strain1}allele.bam.bai \
-g ${strain1}.genome \
-o ${hybrid}_merged_${strain1}allele.HMMRATAC \
--bedgraph TRUE

### Hybrid allele - Strain 2
java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/HMMRATAC_V1.2.5_exe.jar \
-b ${hybrid}_merged_${strain2}allele.bam \
-i ${hybrid}_merged_${strain2}allele.bam.bai \
-g ${strain2}.genome \
-o ${hybrid}_merged_${strain2}allele.HMMRATAC \
--bedgraph TRUE

## Clean up - extract bed columns
awk '{print $1"\t"$7"\t"$8}' ${strain1}_merged_${strain1}allele.HMMRATAC_peaks.gappedPeak \
> ${strain1}_merged_${strain1}allele.HMMRATAC_accessible.bed

awk '{print $1"\t"$7"\t"$8}' ${strain2}_merged_${strain2}allele.HMMRATAC_peaks.gappedPeak \
> ${strain2}_merged_${strain2}allele.HMMRATAC_accessible.bed

awk '{print $1"\t"$7"\t"$8}' ${hybrid}_merged_${strain1}allele.HMMRATAC_peaks.gappedPeak \
> ${hybrid}_merged_${strain1}allele.HMMRATAC_accessible.bed

awk '{print $1"\t"$7"\t"$8}' ${hybrid}_merged_${strain2}allele.HMMRATAC_peaks.gappedPeak \
> ${hybrid}_merged_${strain2}allele.HMMRATAC_accessible.bed

############# Convert strain (1/2) genome -> reference genome 1 -> reference genome 2 ###############
# Strain 1/2 genome/allele -> Reference genome 1
## Strain 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${strain1}_merged_${strain1}allele.HMMRATAC_accessible.bed \
$strain1_to_refgenome1_chain \
${strain1}_merged_${strain1}allele.HMMRATAC_accessible_${ref_genome1}.bed \
${strain1}_merged_${strain1}allele.HMMRATAC_accessible_${ref_genome1}_unmapped

## Strain 2
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${strain2}_merged_${strain2}allele.HMMRATAC_accessible.bed \
$strain2_to_refgenome1_chain \
${strain2}_merged_${strain2}allele.HMMRATAC_accessible_${ref_genome1}.bed \
${strain2}_merged_${strain2}allele.HMMRATAC_accessible_${ref_genome1}_unmapped

## Hybrid allele - Strain 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${hybrid}_merged_${strain1}allele.HMMRATAC_accessible.bed \
$strain1_to_refgenome1_chain \
${hybrid}_merged_${strain1}allele.HMMRATAC_accessible_${ref_genome1}.bed \
${hybrid}_merged_${strain1}allele.HMMRATAC_accessible_${ref_genome1}_unmapped

## Hybrid allele - Strain 2
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${hybrid}_merged_${strain2}allele.HMMRATAC_accessible.bed \
$strain2_to_refgenome1_chain \
${hybrid}_merged_${strain2}allele.HMMRATAC_accessible_${ref_genome1}.bed \
${hybrid}_merged_${strain2}allele.HMMRATAC_accessible_${ref_genome1}_unmapped

# Reference genome 1 -> Reference genome 2
## Strain 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${strain1}_merged_${strain1}allele.HMMRATAC_accessible_dm3.bed \
$refgenome1_to_refgenome2_chain \
${strain1}_merged_${strain1}allele.HMMRATAC_accessible_${ref_genome2}.bed \
${strain1}_merged_${strain1}allele.HMMRATAC_accessible_${ref_genome2}_unmapped

## Strain 2
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${strain2}_merged_${strain2}allele.HMMRATAC_accessible_dm3.bed \
$refgenome1_to_refgenome2_chain \
${strain2}_merged_${strain2}allele.HMMRATAC_accessible_${ref_genome2}.bed \
${strain2}_merged_${strain2}allele.HMMRATAC_accessible_${ref_genome2}_unmapped

## Hybrid allele - Strain 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${hybrid}_merged_${strain1}allele.HMMRATAC_accessible_dm3.bed \
$refgenome1_to_refgenome2_chain \
${hybrid}_merged_${strain1}allele.HMMRATAC_accessible_${ref_genome2}.bed \
${hybrid}_merged_${strain1}allele.HMMRATAC_accessible_${ref_genome2}_unmapped

## Hybrid allele - Strain 2
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${hybrid}_merged_${strain2}allele.HMMRATAC_accessible_dm3.bed \
$refgenome1_to_refgenome2_chain \
${hybrid}_merged_${strain2}allele.HMMRATAC_accessible_${ref_genome2}.bed \
${hybrid}_merged_${strain2}allele.HMMRATAC_accessible_${ref_genome2}_unmapped

################### Make the master peak set and grab reads for each sample #########################
# Concatenate all BED files with common reference genome
cat ${strain1}_merged_${strain1}allele.HMMRATAC_accessible_${ref_genome2}.bed \
${strain2}_merged_${strain2}allele.HMMRATAC_accessible_${ref_genome2}.bed \
${hybrid}_merged_${strain1}allele.HMMRATAC_accessible_${ref_genome2}.bed \
${hybrid}_merged_${strain2}allele.HMMRATAC_accessible_${ref_genome2}.bed \
> ${hybrid}_cat_peaks.bed

# Run mergeBED on concatenated file to create one master peak set
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/./mergeBed -i ${hybrid}_cat_peaks.bed \
> ${hybrid}_merged_peaks.bed

# Get read count for each sample/replicate for master peak set
bedtools multicov -bams ${strain1}_1_${strain1}_sorted_coord.bam \
${strain1}_2_${strain1}_sorted_coord.bam \
${strain1}_3_${strain1}_sorted_coord.bam \
${strain2}_1_${strain2}_sorted_coord.bam \
${strain2}_2_${strain2}_sorted_coord.bam \
${strain2}_3_${strain2}_sorted_coord.bam \
${hybrid}_1_${strain1}_sorted_coord.bam \
${hybrid}_2_${strain1}_sorted_coord.bam \
${hybrid}_3_${strain1}_sorted_coord.bam \
${hybrid}_1_${strain2}_sorted_coord.bam \
${hybrid}_2_${strain2}_sorted_coord.bam \
${hybrid}_3_${strain2}_sorted_coord.bam \
-bed ${hybrid}_merged_peaks_HMMRATAC.bed
> ${hybrid}_merged_HMMRATAC_counts_final.bed
