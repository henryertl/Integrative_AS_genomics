#!/bin/bash

#load module trimgalore

# define folder paths
raw_reads_folder=/nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/RAW_READS
trimmed_reads_and_report_folder=/nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/TRIMMED_READS
fastqc_folder_folder=/nfs/turbo/lsa-wittkopp/Lab/Henry/AS_Integrative_project/FASTQC_REPORT

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
