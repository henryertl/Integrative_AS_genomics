###command to filter reads concordantly mapped only 1x
samtools view -hf 0x2 $1 | grep -v "XS:i:" > $2

##usage ./filter_1x alignments.bam output_filtered_file.sam

###separate chromosomes into two files by allelic origin
samtools view -b in.bam $1 $2 $3 $4 > $5

##usage ./separate_by_chromosomes chr1 chr2 chr3 chr4 zhr_AS_reads.bam
