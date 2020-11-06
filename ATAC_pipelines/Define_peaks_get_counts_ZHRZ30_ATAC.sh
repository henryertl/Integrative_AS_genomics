#!/bin/bash

# define files
strain1=ZHR
strain2=Z30
hybrid=ZHR_Z30
datatype=ATAC
ref_genome1="dm3"
ref_genome1B="sim1"
ref_genome2="dm6"
strain1_to_refgenome1_chain=/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/chainfiles/zhr_bow_to_dm3.chain
strain2_to_refgenome1_chain=/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/chainfiles/z30_bow_to_dm3.chain
refgenome1_to_refgenome2_chain=/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/chainfiles/dm3ToDm6.over.chain
genome2_genome_fasta=/nfs/turbo/lsa-wittkopp/Lab/Henry/genome_fasta_and_index/dm6.fa

# move to working directory
cd /nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/MERGED_AND_FINAL_FILES/${datatype}

################# Define peakset #########################
# Merge replicates for each genotype
## Strain 1

samtools merge -f ${strain1}_merged_${strain1}allele_${datatype}.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_1_${datatype}/${strain1}_1_${datatype}.${strain1}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_2_${datatype}/${strain1}_2_${datatype}.${strain1}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_3_${datatype}/${strain1}_3_${datatype}.${strain1}_sorted_coord.bam

samtools sort -o ${strain1}_merged_${strain1}allele_${datatype}.bam \
${strain1}_merged_${strain1}allele_${datatype}.bam
samtools index ${strain1}_merged_${strain1}allele_${datatype}.bam

## Strain 2
samtools merge -f ${strain2}_merged_${strain2}allele_${datatype}.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_1_${datatype}/${strain2}_1_${datatype}.${strain2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_2_${datatype}/${strain2}_2_${datatype}.${strain2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_3_${datatype}/${strain2}_3_${datatype}.${strain2}_sorted_coord.bam

samtools sort -o ${strain2}_merged_${strain2}allele_${datatype}.bam \
${strain2}_merged_${strain2}allele_${datatype}.bam
samtools index ${strain2}_merged_${strain2}allele_${datatype}.bam

## Hybrid allele - Strain 1
samtools merge -f ${hybrid}_merged_${strain1}allele_${datatype}.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain1}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain1}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain1}_sorted_coord.bam

samtools sort -o ${hybrid}_merged_${strain1}allele_${datatype}.bam \
${hybrid}_merged_${strain1}allele_${datatype}.bam
samtools index ${hybrid}_merged_${strain1}allele_${datatype}.bam

## Hybrid allele - Strain 2
samtools merge -f ${hybrid}_merged_${strain2}allele_${datatype}.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain2}_sorted_coord.bam

samtools sort -o ${hybrid}_merged_${strain2}allele_${datatype}.bam \
${hybrid}_merged_${strain2}allele_${datatype}.bam
samtools index ${hybrid}_merged_${strain2}allele_${datatype}.bam

# Define peaks with HMMRATAC
## Generate genome files needed for HMMRATAC
samtools view -H ${strain1}_merged_${strain1}allele_${datatype}.bam \
| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' \
> ${strain1}.genome

samtools view -H ${strain2}_merged_${strain2}allele_${datatype}.bam \
| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' \
> ${strain2}.genome

## Run HMMRATAC
### Strain 1
java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/HMMRATAC_V1.2.5_exe.jar \
-b ${strain1}_merged_${strain1}allele_${datatype}.bam \
-i ${strain1}_merged_${strain1}allele_${datatype}.bam.bai \
-g ${strain1}.genome \
-o ${strain1}_merged_${strain1}allele_${datatype}.HMMRATAC \
--bedgraph TRUE

### Strain 2
java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/HMMRATAC_V1.2.5_exe.jar \
-b ${strain2}_merged_${strain2}allele_${datatype}.bam \
-i ${strain2}_merged_${strain2}allele_${datatype}.bam.bai \
-g ${strain2}.genome \
-o ${strain2}_merged_${strain2}allele_${datatype}.HMMRATAC \
--bedgraph TRUE

### Hybrid allele - Strain 1
java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/HMMRATAC_V1.2.5_exe.jar \
-b ${hybrid}_merged_${strain1}allele_${datatype}.bam \
-i ${hybrid}_merged_${strain1}allele_${datatype}.bam.bai \
-g ${strain1}.genome \
-o ${hybrid}_merged_${strain1}allele_${datatype}.HMMRATAC \
--bedgraph TRUE

### Hybrid allele - Strain 2
java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/HMMRATAC_V1.2.5_exe.jar \
-b ${hybrid}_merged_${strain2}allele_${datatype}.bam \
-i ${hybrid}_merged_${strain2}allele_${datatype}.bam.bai \
-g ${strain2}.genome \
-o ${hybrid}_merged_${strain2}allele_${datatype}.HMMRATAC \
--bedgraph TRUE

## Clean up - extract bed columns
awk '{print $1"\t"$7"\t"$8}' ${strain1}_merged_${strain1}allele_${datatype}.HMMRATAC_peaks.gappedPeak \
> ${strain1}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible.bed

awk '{print $1"\t"$7"\t"$8}' ${strain2}_merged_${strain2}allele_${datatype}.HMMRATAC_peaks.gappedPeak \
> ${strain2}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible.bed

awk '{print $1"\t"$7"\t"$8}' ${hybrid}_merged_${strain1}allele_${datatype}.HMMRATAC_peaks.gappedPeak \
> ${hybrid}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible.bed

awk '{print $1"\t"$7"\t"$8}' ${hybrid}_merged_${strain2}allele_${datatype}.HMMRATAC_peaks.gappedPeak \
> ${hybrid}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible.bed

######### Convert strain (1/2) genome -> reference genome 1 -> reference genome 2 #########
## Strain 1 genome-> Reference genome 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${strain1}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible.bed \
$strain1_to_refgenome1_chain \
${strain1}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_${ref_genome1}.bed \
${strain1}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_${ref_genome1}_unmapped

## Strain 2 -> Ref genome 1B
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${strain2}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible.bed \
$strain2_to_refgenome1B_chain \
${strain2}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_${ref_genome1B}.bed \
${strain2}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_${ref_genome1B}_unmapped

## Hybrid allele - Strain 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${hybrid}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible.bed \
$strain1_to_refgenome1_chain \
${hybrid}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_${ref_genome1}.bed \
${hybrid}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_${ref_genome1}_unmapped

## Hybrid allele - Strain 2
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${hybrid}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible.bed \
$strain2_to_refgenome1_chain \
${hybrid}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_${ref_genome1}.bed \
${hybrid}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_${ref_genome1}_unmapped

# Reference genome 1 -> Reference genome 2
## Strain 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${strain1}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_dm3.bed \
$refgenome1_to_refgenome2_chain \
${strain1}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}.bed \
${strain1}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}_unmapped

## Strain 2
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${strain2}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_dm3.bed \
$refgenome1_to_refgenome2_chain \
${strain2}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}.bed \
${strain2}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}_unmapped

## Hybrid allele - Strain 1
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${hybrid}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_dm3.bed \
$refgenome1_to_refgenome2_chain \
${hybrid}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}.bed \
${hybrid}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}_unmapped

## Hybrid allele - Strain 2
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${hybrid}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_dm3.bed \
$refgenome1_to_refgenome2_chain \
${hybrid}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}.bed \
${hybrid}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}_unmapped

################### Make the master peak set #########################
# Concatenate all BED files with common reference genome
cat ${strain1}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}.bed \
${strain2}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}.bed \
${hybrid}_merged_${strain1}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}.bed \
${hybrid}_merged_${strain2}allele_${datatype}.HMMRATAC_accessible_${ref_genome2}.bed \
> ${hybrid}_${datatype}_cat_peaks.bed

# Run mergeBED on concatenated file to create one master peak set
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/./mergeBed -i ${hybrid}_${datatype}_cat_peaks.bed \
> ${hybrid}_${datatype}_merged_peaks.bed

# Get read count for each sample/replicate for master peak set
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/./multiBamCov -bams ../../ALL_${datatype}_ALIGNED/${strain1}_1_${datatype}/${strain1}_1_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_2_${datatype}/${strain1}_2_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain1}_3_${datatype}/${strain1}_3_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_1_${datatype}/${strain2}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_2_${datatype}/${strain2}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${strain2}_3_${datatype}/${strain2}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain1}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_1_${datatype}/${hybrid}_1_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_2_${datatype}/${hybrid}_2_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
../../ALL_${datatype}_ALIGNED/${hybrid}_3_${datatype}/${hybrid}_3_${datatype}.${strain2}_${ref_genome2}_sorted_coord.bam \
-bed ${hybrid}_${datatype}_merged_peaks.bed \
> ${hybrid}_${datatype}_merged_counts_final.bed



java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/HMMRATAC_V1.2.5_exe.jar \
-b Z30_merged_Z30allele_ATAC.bam \
-i Z30_merged_Z30allele_ATAC.bam.bai \
-g Z30.genome \
-o Z30_merged_Z30allele_ATAC.HMMRATAC \
--bedgraph TRUE \
--upper 40
