#!/bin/bash

# how to use script
## if on cluster, load modules: samtools, bowtie2
## ./ZHR_Z30_align_reads_and_convert_to_ref.sh [desired file prefix] /scratch/lsa_root/lsa/hertl/genomes/indexed_genomes/zhr_z30_bow_cat_GATC_indexed /path/to/forward_reads /path/to/reverse_reads

# define liftover chain files and genome file for conversions
zhr_to_dm3_chain=/scratch/lsa_root/lsa/hertl/liftover_utilities/tempo_mode_genomes_chainfiles/Dmel_zhr/zhr_bow_to_dm3.chain
z30_to_dm3_chain=/scratch/lsa_root/lsa/hertl/liftover_utilities/tempo_mode_genomes_chainfiles/Dmel_z30/z30_bow_to_dm3.chain
dm3_to_dm6_chain=/scratch/lsa_root/lsa/hertl/other_datasets/dros_sp_EA_ATAC_raw/liftover_chain_files/dm3ToDm6.over.chain
dm6_genome_fasta=/scratch/lsa_root/lsa/hertl/genomes/dm6.fa
indexed_genome=/path/to/indexed/genome
forward_reads=/path/to/reads
reverse_reads=/path/to/reads
SAMPLE_NAME={sample_prefix}_atac

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

# convert BAM to fastq pe reads to re-align to cat indexed_genome
## output gives {SAMPLE_NAME}_1_WASP_filtered.fq.gz and {SAMPLE_NAME}_2_WASP_filtered.fq.gz
bedtools bamtofastq -i merge/${SAMPLE_NAME}.keep.merge.sort_name.rmdups.bam \
        -fq /scratch/lsa_root/lsa/hertl/AS_ATAC_new/${SAMPLE_NAME}_1_WASP_filtered.fq.gz \
        -fq2 /scratch/lsa_root/lsa/hertl/AS_ATAC_new/[${SAMPLE_NAME}_2_WASP_filtered.fq.gz

######################################################################
#### align WASP filtered reads to concatenated genome and AS sort ####
######################################################################

# change directory to strain sp
cd /scratch/lsa_root/lsa/hertl/AS_ATAC_new/[strain_prefix]

############# align reads, clean up, and sort alleles ################
# align reads
bowtie2 --very-sensitive -x $indexed_genome \
        -1 ${SAMPLE_NAME}_1_WASP_filtered.fq.gz \
        -2 ${SAMPLE_NAME}_2_WASP_filtered.fq.gz \
        -S ${SAMPLE_NAME}.sam
samtools sort -O BAM -o ${SAMPLE_NAME}.bam \
        ${SAMPLE_NAME}.sam
samtools index ${SAMPLE_NAME}.bam

##parent 1
samtools view -b -h -f 3 -F 4 -F 256 -F 1024 -F 2048 -q 20 \
        ${SAMPLE_NAME}.bam zhr_bow_chr2L zhr_bow_chr2R zhr_bow_chr3L zhr_bow_chr3R zhr_bow_chr4 zhr_bow_chrX > ${SAMPLE_NAME}.ZHR_sorted.bam
##parent 2
samtools view -b -h -f 3 -F 4 -F 256 -F 1024 -F 2048 -q 20 \
        ${SAMPLE_NAME}.bam z30_bow_chr2L z30_bow_chr2R z30_bow_chr3L z30_bow_chr3R z30_bow_chr4 z30_bow_chrX > ${SAMPLE_NAME}.Z30_sorted.bam

################### convert to reference genomes ####################

# sort by name and coordinate
samtools sort -n -o ${SAMPLE_NAME}.ZHR_sorted_name.bam \
        ${SAMPLE_NAME}.ZHR_sorted.bam
samtools sort -n -o ${SAMPLE_NAME}.Z30_sorted_name.bam \
        ${SAMPLE_NAME}.Z30_sorted.bam

# convert bam to bed with bedtools
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/bamToBed -cigar -i ${SAMPLE_NAME}.ZHR_sorted_name.bam > ${SAMPLE_NAME}.ZHR_sorted_name.bed
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/bamToBed -cigar -i ${SAMPLE_NAME}.Z30_sorted_name.bam > ${SAMPLE_NAME}.Z30_sorted_name.bed

# convert from strain sp to dm3
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${SAMPLE_NAME}.ZHR_sorted_name.bed $zhr_to_dm3_chain ${SAMPLE_NAME}.ZHR_dm3.bed ${SAMPLE_NAME}.ZHR_dm3_unmapped
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${SAMPLE_NAME}.Z30_sorted_name.bed $z30_to_dm3_chain ${SAMPLE_NAME}.Z30_dm3.bed ${SAMPLE_NAME}.Z30_dm3_unmapped

# convert dm3 to dm6
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${SAMPLE_NAME}.ZHR_dm3.bed $dm3_to_dm6_chain ${SAMPLE_NAME}.ZHR_dm6.bed ${SAMPLE_NAME}.ZHR_dm6_unmapped
/nfs/turbo/lsa-wittkopp/Lab/Henry/liftover_utilities/./liftOver.dms ${SAMPLE_NAME}.Z30_dm3.bed $dm3_to_dm6_chain ${SAMPLE_NAME}.Z30_dm6.bed ${SAMPLE_NAME}.Z30_dm6_unmapped

# convert dm6 beds to dm6 bam
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/bedToBam -i ${SAMPLE_NAME}.ZHR_dm6.bed -g $dm6_genome_fasta > ${SAMPLE_NAME}.ZHR_dm6.bam
/nfs/turbo/lsa-wittkopp/Lab/Henry/executables/bedtools2/bin/bedToBam -i ${SAMPLE_NAME}.Z30_dm6.bed -g $dm6_genome_fasta > ${SAMPLE_NAME}.Z30_dm6.bam

# sort and index final file
samtools sort -o ${SAMPLE_NAME}.ZHR_dm6_sorted.bam \
        ${SAMPLE_NAME}.ZHR_dm6.bam
samtools index ${SAMPLE_NAME}.ZHR_dm6_sorted.bam

samtools sort -o ${SAMPLE_NAME}.Z30_dm6_sorted.bam \
        ${SAMPLE_NAME}.Z30_dm6.bam
samtools index ${SAMPLE_NAME}.Z30_dm6_sorted.bam

# sort originals by coordinate for next step and move this file to merged_files directory for next step
samtools sort -O BAM -o ${SAMPLE_NAME}.ZHR_sorted_coord.bam \
        ${SAMPLE_NAME}.ZHR_sorted.bam
mv ${SAMPLE_NAME}.ZHR_sorted_coord.bam /scratch/lsa_root/lsa/hertl/AS_ATAC_new/ZHR_Z30_gene_counts

samtools sort -O BAM -o ${SAMPLE_NAME}.Z30_sorted_coord.bam \
        ${SAMPLE_NAME}e.Z30_sorted.bam
mv ${SAMPLE_NAME}.Z30_sorted_coord.bam /scratch/lsa_root/lsa/hertl/AS_ATAC_new/ZHR_Z30_gene_counts
