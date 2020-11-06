#!/bin/bash

# define files/
indexed_genome=/scratch/lsa_root/lsa/hertl/genomes/indexed_genomes/zhr_bow_GATC_z30_snp_masked
forward_reads=/scratch/lsa_root/lsa/hertl/other_datasets/tempo_mode_gDNA/Z30_ZHR_hyb_peRNAseq_F.fq
reverse_reads=/scratch/lsa_root/lsa/hertl/other_datasets/tempo_mode_gDNA/Z30_ZHR_hyb_peRNAseq_R.fq
SAMPLE_NAME=zhr_z30_temponmode

############# WASP #################
## assuming have already generated SNP bed file (see step 1 from WASP github), and SNP files, and have indexed masked genome
## make sure forward and reverse reads are named {SAMPLE_NAME_1/2.fq.gz}
####################################

cd /scratch/lsa_root/lsa/hertl/executables/WASP

# align to masked genome
bowtie2 -x $indexed_genome_masked -1 $forward_reads \
           -2 $reverse_reads \
           | samtools view -b -q 10 - > map1/${SAMPLE_NAME}.bam
   samtools sort -o map1/${SAMPLE_NAME}.sort.bam map1/${SAMPLE_NAME}.bam
   samtools index map1/${SAMPLE_NAME}.sort.bam

# find reads from above BAM file intersecting snps
python mapping/find_intersecting_snps.py \
         --is_paired_end \
         --is_sorted \
         --output_dir /find_intersecting_snps \
         --snp_dir /scratch/lsa_root/lsa/hertl/genomes/SNP_files/WASP_zhr_z30 \
         map1/${SAMPLE_NAME}.sort.bam

# remap reads that overlap SNPs (the remap reads output from above)
bowtie2 -x $indexed_genome_masked \
             -1 find_intersecting_snps/${SAMPLE_NAME}_1.remap.fq.gz \
             -2 find_intersecting_snps/${SAMPLE_NAME}_2.remap.fq.gz \
         | samtools view -b -q 10 - > map2/${SAMPLE_NAME}.bam
   samtools sort -o map2/${SAMPLE_NAME}.sort.bam map2/${SAMPLE_NAME}.bam
   samtools index map2/${SAMPLE_NAME}.sort.bam

# filter to remove reads that don't remap same location after allele swap
python mapping/filter_remapped_reads.py \
     find_intersection_snps/${SAMPLE_NAME}.to.remap.bam \
     map2/${SAMPLE_NAME}.sort.bam \
     filter_remapped_reads/${SAMPLE_NAME}.keep.bam

# merge things
samtools merge merge/${SAMPLE_NAME}.keep.merge.bam \
         filter_remapped_reads/${SAMPLE_NAME}.keep.bam  \
         find_intersecting_snps/${SAMPLE_NAME}.keep.bam
samtools sort -o  merge/${SAMPLE_NAME}.keep.merge.sort.bam \
         merge/${SAMPLE_NAME}.keep.merge.bam
samtools index ${SAMPLE_NAME}.keep.merged.sort.bam

# remove duplicates from final BAM
python rmdup_pe.py merge/${SAMPLE_NAME}.keep.merge.sort.bam merge/${SAMPLE_NAME}.keep.merge.sort.rmdups.bam
samtools sort -o merge/${SAMPLE_NAME}.keep.merge.sort_name.rmdups.bam merge/${SAMPLE_NAME}.keep.merge.sort.rmdups.bam
samtools index merge/${SAMPLE_NAME}.keep.merge.sort_name.rmdups.bam
