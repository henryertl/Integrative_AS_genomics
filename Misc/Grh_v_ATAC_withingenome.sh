### GOAL: correlation between Grh binding intensity and ATAC signal/accessiblity?

# ./ATAC_v_bed [Grh_input] [ATAC_input] [output_prefix]

# Generate file with merged Grh binding and ATAC signal
## Run bedtools intersect


bedtools intersect -wa -a $1 -b $2 > $3.bed

Rscript /path/to/Rscript

bedtools intersect -wao -a /Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ZHR_merged_ZHRallele_ATAC.HMMRATAC_peaks_dm3_test.gappedPeak -b /Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/Grh_data/GSM2716905_ChIPmentation_grh-GFP_anti-GFP_Wing.bed  > Grh_ATAC_ZHR_Z30.bed
