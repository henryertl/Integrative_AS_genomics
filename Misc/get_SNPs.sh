#!/bin/bash

#################### Goal: get SNP (vcf and bed) files from raw reads #########################
########### This script assumes you have raw gDNA reads with adapters trimmed off #############
############################ and assembled genomes for both strains ###########################

# Define inputs
strain1=ZHR
strain2=Z30
strain_1_genome_path=/nfs/turbo/lsa-wittkopp/Lab/Henry/genome_fasta_and_index/ZHR_genome.fa
strain_2_genome_path=/nfs/turbo/lsa-wittkopp/Lab/Henry/genome_fasta_and_index/Z30_genome.fa

# Change to working directory with reads
cd /nfs/turbo/lsa-wittkopp/Lab/Henry/other_datasets/tempo_mode_gDNA

# Index genomes (ignore if already done)
#bowtie2-build $strain_1_genome_path ${strain1}_genome_indexed
#bowtie2-build $strain_2_genome_path ${strain2}_genome_indexed

######## align reads to genome ##########
## Strain 2 to Strain 3
bowtie2 --very-sensitive -x ${strain2}_genome_indexed \
        -1 ${strain1}_gDNA_F_cutadapt.fastq.gz \
        -2 ${strain1}_gDNA_R_cutadapt.fastq.gz \
        -S ${strain1}_gDNA_aligned_to_${strain2}.sam

samtools sort -O BAM -o ${strain1}_gDNA_aligned_to_${strain2}.bam \
        ${strain1}_gDNA_aligned_to_${strain2}.sam
samtools index ${strain1}_gDNA_aligned_to_${strain2}.bam

## Strain 2 to Strain 1
bowtie2 --very-sensitive -x ${strain1}_genome_indexed \
        -1 ${strain2}_gDNA_F_cutadapt.fastq.gz \
        -2 ${strain2}_gDNA_R_cutadapt.fastq.gz \
        -S ${strain2}_gDNA_aligned_to_${strain1}.sam

samtools sort -O BAM -o ${strain2}_gDNA_aligned_to_${strain1}.bam \
        ${strain2}_gDNA_aligned_to_${strain1}.sam
samtools index ${strain2}_gDNA_aligned_to_${strain1}.bam

######### Call SNPs with GATK Haplotype caller #########

# GATK requires BAM files to have read groups (which they do not), so have to artifically add in
## Strain 1 to Strain 2
java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/picard-2.jar AddOrReplaceReadGroups \
I=${strain1}_gDNA_aligned_to_${strain2}.bam \
O=${strain1}_gDNA_aligned_to_${strain2}_RGadded.bam \
RGID=4 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=20

### Sort, index, and replace header
samtools sort -O BAM -o ${strain1}_gDNA_aligned_to_${strain2}_RGadded.bam ${strain1}_gDNA_aligned_to_${strain2}_RGadded.bam
samtools index ${strain1}_gDNA_aligned_to_${strain2}_RGadded.bam
samtools view -H ${strain1}_gDNA_aligned_to_${strain2}.bam |
  samtools reheader - ${strain1}_gDNA_aligned_to_${strain2}_RGadded.bam

## Strain 2 to Strain 1
java -jar /nfs/turbo/lsa-wittkopp/Lab/Henry/executables/picard-2.jar AddOrReplaceReadGroups \
I=${strain2}_gDNA_aligned_to_${strain1}.bam \
O=${strain2}_gDNA_aligned_to_${strain1}_RGadded.bam \
RGID=4 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=20

### sort, index, and replace header
samtools sort -O BAM -o ${strain2}_gDNA_aligned_to_${strain1}_RGadded.bam ${strain2}_gDNA_aligned_to_${strain1}_RGadded.bam
samtools index ${strain2}_gDNA_aligned_to_${strain1}_RGadded.bam
samtools view -H ${strain2}_gDNA_aligned_to_${strain1}.bam |
  samtools reheader - ${strain2}_gDNA_aligned_to_${strain1}_RGadded.bam
# Run GATK HaplotypeCaller to generate vcf file

## Strain 1 to Strain 2
gatk HaplotypeCaller -I ${strain1}_gDNA_aligned_to_${strain2}_RGadded.bam \
  -R $strain_2_genome_path \
  -O ${strain1}_gDNA_aligned_to_${strain2}_SNPs.vcf \
  --genotyping-mode DISCOVERY \
  --output-mode EMIT_ALL_SITES \
  --standard-min-confidence-threshold-for-calling 30

## Strain 2 to Strain 1
gatk HaplotypeCaller -I ${strain2}_gDNA_aligned_to_${strain1}_RGadded.bam \
  -R $strain_1_genome_path \
  -O ${strain2}_gDNA_aligned_to_${strain1}_SNPs.vcf \
  --genotyping-mode DISCOVERY \
  --output-mode EMIT_ALL_SITES \
  --standard-min-confidence-threshold-for-calling 30

# Use BEDops software to convert VCF to BED for downstream analyses

## Strain 1 to Strain 2
vcf2bed --deletions < ${strain1}_gDNA_aligned_to_${strain2}_SNPs.vcf > ${strain1}_gDNA_aligned_to_${strain2}_deletions.bed
vcf2bed --snvs < ${strain1}_gDNA_aligned_to_${strain2}_SNPs.vcf > ${strain1}_gDNA_aligned_to_${strain2}_SNPs_snvs.bed
bedops --everything ${strain1}_gDNA_aligned_to_${strain2}_SNPs_{deletions,snvs}.bed > ${strain1}_gDNA_aligned_to_${strain2}_snvs_deletions.bed

## Strain 2 to Strain 1
vcf2bed --deletions < ${strain2}_gDNA_aligned_to_${strain1}_SNPs.vcf > ${strain2}_gDNA_aligned_to_${strain1}_deletions.bed
vcf2bed --snvs < ${strain2}_gDNA_aligned_to_${strain1}_SNPs.vcf > ${strain2}_gDNA_aligned_to_${strain1}_SNPs_snvs.bed
bedops --everything ${strain2}_gDNA_aligned_to_${strain1}_SNPs_{deletions,snvs}.bed > ${strain2}_gDNA_aligned_to_${strain1}_snvs_deletions.bed

# Use BEDtools to compile both directions (strain 1 - strain 2)


# Clean up and organinze files
