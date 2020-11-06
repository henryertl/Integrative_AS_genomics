# convert bam to bed then get counts from sliding window
##D_wil
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_1_dm3_md_hqaa_sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_2_dm3_md_hqaa_sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_wil/D_wil_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

##D_mel
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_1_dm3_md_hqaa_sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_2_dm3_md_hqaa_sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_mel/D_mel_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

##D_sim
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_1_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_2_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sim/D_sim_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

##D_per
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_1_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_2_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_per/D_per_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

##D_vir
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_1_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_2_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_vir/D_vir_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

##D_yak
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_1_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_2_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_yak/D_yak_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

##D_pse
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_1_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_2_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_pse/D_pse_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

##D_ere
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_1_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_2_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ere/D_ere_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

##D_moj
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_1_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_2_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_moj/D_moj_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

##D_sech
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_1_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_2_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_sech/D_sech_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

##D_ana
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_1_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_1_dm3_md_hqaa_sorted.bed
convert2bed --input=bam [--output=bed] < /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_2_dm3_md_hqaa.sorted.bam > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_2_dm3_md_hqaa_sorted.bed

bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_1_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed
bedmap --faster --echo --count /scratch/lsa_flux/hertl/executables/bedops_bin/dm3_windows_1000.bed /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_2_dm3_md_hqaa_sorted.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed

sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_1_dm3_md_hqaa_sorted_1000_windows_final.bed
sed 's/|/\t/g' /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows.bed > /scratch/lsa_flux/hertl/dros_sp_EA_ATAC_raw/alignment_files/D_ana/D_ana_EA_ATAC_2_dm3_md_hqaa_sorted_1000_windows_final.bed

## merge all final files
