#!/bin/bash

#load module trimgalore

# define folder paths
raw_reads_folder=/scratch/lsa_root/lsa/hertl/1777-HE/all_raw_reads/ATAC_raw
trimmed_reads_and_report_folder=/scratch/lsa_root/lsa/hertl/1777-HE/TRIMMED_READS
fastqc_folder_folder=/scratch/lsa_root/lsa/hertl/1777-HE/FASTQC_REPORT

# cd to folder with raw reads
cd $raw_reads_folder

# trim all Nextera (replace if necessary)
parallel --xapply \
trim_galore \
--nextera \
--paired \
--fastqc \
-o ./ \
::: *_R1_001.fastq.gz ::: *_R2_001.fastq.gz

# clean up and organize files
mv *_val_* $trimmed_reads_and_report_folder
mv *fastqc* $fastqc_folder_folder
