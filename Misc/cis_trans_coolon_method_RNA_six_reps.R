###############################################################################
################# Coolon Method for ASE cis - trans analysis ##################
###############################################################################

######################### hypergeometric sampling #############################

# read in dataset
df <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/ZHR_Z30_genic_counts_combined_final_dm6.txt")

# get total read counts
sum(df[,2])
# 24837936
sum(df[,3])
# 8142640
sum(df[,4])
# 6555516
sum(df[,5])
# 2297297
sum(df[,6])
# 975085
sum(df[,7])
# 278809
sum(df[,8])
# 653762
sum(df[,9])
# 168137
sum(df[,10])
# 24837936
sum(df[,11])
# 8142640
sum(df[,12])
# 6555516
sum(df[,13])
# 2297297
sum(df[,14])
# 975085
sum(df[,15])
# 278809
sum(df[,16])
# 653762
sum(df[,17])
# 168137
sum(df[,18])
# 653762
min <- sum(df[,19])
# 168137
sum(df[,20])
# 24837936
sum(df[,21])
# 8142640
sum(df[,22])
# 6555516
sum(df[,23])
# 2297297
sum(df[,24])
# 975085
sum(df[,25])
# 278809

# randomly down sample each to 227904 reads and add to new file
df$number <- (1:nrow(df))

df_downsamp <- NULL %>% as.data.frame()
df_downsamp <- df[c(1:3)]
df_downsamp$ID <- cbind(1:nrow(df)) %>% as.data.frame()

######### P1_1 ############
# expand read counts and sample by min number
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,2])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P1_1_downsamp")
d_merge <- merge(d_table, df, by.x = "ID", by.y = "number")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)

######### P1_2 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,3])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P1_2_downsamp")
d_merge2 <- merge(d_table, d_merge, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge)

######### P1_3 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,4])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P1_3_downsamp")
d_merge3 <- merge(d_table, d_merge2, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge2)

######### P1_4 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,5])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P1_4_downsamp")
d_merge4 <- merge(d_table, d_merge3, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge3)

######### P1_5 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,6])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P1_5_downsamp")
d_merge5 <- merge(d_table, d_merge4, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge4)

######### P1_6 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,7])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P1_6_downsamp")
d_merge6 <- merge(d_table, d_merge5, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge5)

######### P2_1 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,8])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P2_1_downsamp")
d_merge7 <- merge(d_table, d_merge6, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge6)

######### P2_2 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,9])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P2_2_downsamp")
d_merge8 <- merge(d_table, d_merge7, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge7)

######### P2_3 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,10])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P2_3_downsamp")
d_merge9 <- merge(d_table, d_merge8, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge8)

######### P2_4 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,11])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P2_4_downsamp")
d_merge10 <- merge(d_table, d_merge9, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge9)

######### P2_5 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,12])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P2_5_downsamp")
d_merge11 <- merge(d_table, d_merge10, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge10)

######### P2_6 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,13])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "P2_6_downsamp")
d_merge12 <- merge(d_table, d_merge11, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge11)

######### HYB_1_P1 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,14])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_1_P1_downsamp")
d_merge13 <- merge(d_table, d_merge12, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge12)

######### HYB_2_P1 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,15])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_2_P1_downsamp")
d_merge14 <- merge(d_table, d_merge13, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge13)

######### HYB_3_P1 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,16])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_3_P1_downsamp")
d_merge15 <- merge(d_table, d_merge14, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge14)

######### HYB_4_P1 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,17])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_4_P1_downsamp")
d_merge16 <- merge(d_table, d_merge15, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge15)

######### HYB_5_P1 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,18])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_5_P1_downsamp")
d_merge17 <- merge(d_table, d_merge16, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge16)

######### HYB_6_P1 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,19])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_6_P1_downsamp")
d_merge18 <- merge(d_table, d_merge17, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge17)

######### HYB_1_P2 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,20])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_1_P2_downsamp")
d_merge19 <- merge(d_table, d_merge18, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge18)

######### HYB_2_P2 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,21])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_2_P2_downsamp")
d_merge20 <- merge(d_table, d_merge19, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge19)

######### HYB_3_P2 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,22])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_3_P2_downsamp")
d_merge21 <- merge(d_table, d_merge20, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge20)

######### HYB_4_P2 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,23])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_4_P2_downsamp")
d_merge22 <- merge(d_table, d_merge21, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge21)

######### HYB_5_P2 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,24])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_5_P2_downsamp")
d_merge23 <- merge(d_table, d_merge22, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge22)

######### HYB_6_P2 ############
a1<-c()
for (i in df[,26]) {
  temp<-rep.int(i,df[i,25])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_6_P2_downsamp")
d_merge24 <- merge(d_table, d_merge23, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge23)

d_merge24 <- as.data.frame(d_merge24)
down_samp_final <- d_merge24[rev(2:26)]
colnames(down_samp_final) <- c("gene", "P1_1", "P1_2", "P1_3", "P1_4", "P1_5", "P1_6",
"P2_1", "P2_2", "P2_3", "P2_4", "P2_5", "P2_6",
"HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_4_P1", "HYB_5_P1", "HYB_6_P1",
"HYB_1_P2", "HYB_2_P2", "HYB_3_P2", "HYB_4_P2", "HYB_5_P2", "HYB_6_P2")
write.table(down_samp_final, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/ZHR_Z30_genic_counts_combined_final_dm6_downsampled.txt", sep = "\t", row.names = F, quote = F)


################ Remove any rows with reads <20 and get means across reps ##################
## this is going to yield different data than the model, which subsets <= 20 with no down sampling ##

ZHR_Z30_reps_merged_dm3_final_downsamp_great20 <- subset(d_merge8, d_merge8$ZHR_1_downsamp >= 20 & d_merge8$ZHR_3_downsamp >= 20 &  d_merge8$Z30_1_downsamp >= 20 &  d_merge8$Z30_2_downsamp >= 20 & d_merge8$HYB_1_ZHR_downsamp >= 20 & d_merge8$HYB_1_Z30_downsamp >= 20 & d_merge8$HYB_2_ZHR_downsamp >= 20 & d_merge8$HYB_2_Z30_downsamp >= 20)

############## make file with only downsampled columns for separate analysis ################

ZHR_Z30_reps_downsamp_great20_final <- NULL
ZHR_Z30_reps_downsamp_great20_final <- cbind(ZHR_Z30_reps_merged_dm3_final_downsamp_great20[10:12], ZHR_Z30_reps_merged_dm3_final_downsamp_great20[2:9])
write.table(ZHR_Z30_reps_downsamp_great20_final, file = "ZHR_Z30_reps_downsamp_great20_final.txt", sep = "\t", row.names = F, col.names = T)

##### means ######

ZHR_Z30_reps_merged_dm3_final_downsamp_great20$ZHR_merged=rowMeans(ZHR_Z30_reps_merged_dm3_final_downsamp_great20[,c("ZHR_1_downsamp", "ZHR_3_downsamp")], na.rm=TRUE)
ZHR_Z30_reps_merged_dm3_final_downsamp_great20$Z30_merged=rowMeans(ZHR_Z30_reps_merged_dm3_final_downsamp_great20[,c("Z30_1_downsamp", "Z30_2_downsamp")], na.rm=TRUE)
ZHR_Z30_reps_merged_dm3_final_downsamp_great20$HYB_ZHR_merged=rowMeans(ZHR_Z30_reps_merged_dm3_final_downsamp_great20[,c("HYB_1_ZHR_downsamp", "HYB_2_ZHR_downsamp")], na.rm=TRUE)
ZHR_Z30_reps_merged_dm3_final_downsamp_great20$HYB_Z30_merged=rowMeans(ZHR_Z30_reps_merged_dm3_final_downsamp_great20[,c("HYB_1_Z30_downsamp", "HYB_2_Z30_downsamp")], na.rm=TRUE)

## write table to be in same format as data2 file from coolon paper
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final <- NULL
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final$chrom <- ZHR_Z30_reps_merged_dm3_final_downsamp_great20$chrom
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final$start <- ZHR_Z30_reps_merged_dm3_final_downsamp_great20$start
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final$end <- ZHR_Z30_reps_merged_dm3_final_downsamp_great20$end
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final$ZHR_merged <- ZHR_Z30_reps_merged_dm3_final_downsamp_great20$ZHR_merged
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final$Z30_merged <- ZHR_Z30_reps_merged_dm3_final_downsamp_great20$Z30_merged
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final$HYB_ZHR_merged <- ZHR_Z30_reps_merged_dm3_final_downsamp_great20$HYB_ZHR_merged
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final$HYB_Z30_merged <- ZHR_Z30_reps_merged_dm3_final_downsamp_great20$HYB_Z30_merged

ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final <- as.data.frame(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final)
write.table(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final, file = "ZHR_Z30_reps_merged_dm3_final_downsamp_great20_final.txt")

################ cis trans stats #################
## round to nearest integer for binomial tests
data2 <- NULL
data.output5 <- NULL
data2$chrom <- ZHR_Z30_reps_merged_dm3_final_downsamp_great20$chrom
data2$start <- ZHR_Z30_reps_merged_dm3_final_downsamp_great20$start
data2$end <- ZHR_Z30_reps_merged_dm3_final_downsamp_great20$end
data2$ZHR_merged <- round(ZHR_Z30_reps_merged_dm3_final_downsamp_great20$ZHR_merged)
data2$Z30_merged <- round(ZHR_Z30_reps_merged_dm3_final_downsamp_great20$Z30_merged)
data2$HYB_ZHR_merged <- round(ZHR_Z30_reps_merged_dm3_final_downsamp_great20$HYB_ZHR_merged)
data2$HYB_Z30_merged <- round(ZHR_Z30_reps_merged_dm3_final_downsamp_great20$HYB_Z30_merged)
data2 <- as.data.frame(data2)

## load packages
library(Hmisc)

## initiate variables
# initiate variables (some won't be used)
pvalsHyb1 <- NULL;
pvalsHyb1_Par <- NULL;
pvalsHyb2 <- NULL;
pvalsHyb2_Par <- NULL;
pvalsPar1 <- NULL;
pvalsPar2 <- NULL;
pvalsHyb1_Hyb2 <- NULL;
pvalsHyb1_Hyb2_FET <- NULL;
pvalsHyb1_Hyb2_FET2 <- NULL;
pvalsHyb1_Bi <- NULL;
pvalsHyb2_Bi <- NULL;

Par1_FET2 <- NULL;
Par2_FET2 <- NULL;

tempH1_H2_2 <- NULL;
tempH1_all <- NULL;
tempH2_all <- NULL;
tempH1 <- NULL;
tempH1w <- NULL;
tempH1_P <- NULL;
tempH2 <- NULL;
tempH2w <- NULL;
tempH2_P <- NULL;
tempP1 <- NULL;
tempP1w <- NULL;
tempH1_H2 <- NULL;
tempP2 <- NULL;
tempP2w <- NULL;
tempH1_H2_all <- NULL;
tempH1H2w <-NULL

for (j in 1:nrow(data2))

{
  Par2_Bi <- binom.test(data2[[j,4]], (data2[[j,4]]+data2[[j,5]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
  #Par2_FET <- fisher.test(matrix(c(data2[[j,2]],(sumC-data2[[j,2]]),data2[[j,3]],(sumD-data2[[j,3]])),nr=2));

  Hyb2_Bi <- binom.test(data2[[j,6]], (data2[[j,6]]+data2[[j,7]]), p = 0.5, alternative = c("t"), conf.level = 0.95);

  # collect p-values from binomial tests
  pvalsHyb2 <- rbind(pvalsHyb2,Hyb2_Bi$p.value);
  pvalsPar2 <- rbind(pvalsPar2,Par2_Bi$p.value);
  #pvalsPar2 <- rbind(pvalsPar2, Par2_FET$p.value);


  tempH2 <- rbind(tempH2,c(Hyb2_Bi$estimate,Hyb2_Bi$conf.int[1],Hyb2_Bi$conf.int[2],Hyb2_Bi$p.value));
  tempP2 <- rbind(tempP2,c(Par2_Bi$estimate,Par2_Bi$conf.int[1],Par2_Bi$conf.int[2],Par2_Bi$p.value));
  #tempP2 <- rbind(tempP2,c(Par2_FET$estimate,Par2_FET$conf.int[1],Par2_FET$conf.int[2],Par2_FET$p.value));

  # collect Wilson confidence intervals for binomial tests (to be used later for plotting seq. vs pyro)
  Par2_Bi2 <- binconf(data2[[j,4]], (data2[[j,4]]+data2[[j,5]]), alpha = 0.05, method = c("wilson"));

  #create space filling columns (#11 and #12) to keep spacing constant for mel-mel and other comparisons(mel-sim, sim-sec), bascially a repeat of 7 and 8
  #Par2_FET2 <-rbind(Par2_FET2,c(Par2_FET$conf.int[1],Par2_FET$conf.int[2]));

  Hyb2_Bi2 <- binconf(data2[[j,6]], (data2[[j,6]]+data2[[j,7]]), alpha = 0.05, method = c("wilson"));

  tempH2w <- rbind(tempH2w,c(Hyb2_Bi2[2],Hyb2_Bi2[3]));
  tempP2w <- rbind(tempP2w,c(Par2_Bi2[2],Par2_Bi2[3]));

  # Fisher's exact test from 2x2 tables of counts
  hyb2_par.FET <- fisher.test(matrix(c(data2[[j,6]],data2[[j,7]],data2[[j,4]],data2[[j,5]]),nr=2));

  # collect p-values from FETs
  pvalsHyb2_Par <- rbind(pvalsHyb2_Par,hyb2_par.FET$p.value);

  tempH2_P <- rbind(tempH2_P,c(hyb2_par.FET$estimate,hyb2_par.FET$conf.int[1],hyb2_par.FET$conf.int[2],hyb2_par.FET$p.value));
}

#FDR correct pvalues
pvalsHyb2.adj <- p.adjust(pvalsHyb2,method="BY");
pvalsPar2.adj <- p.adjust(pvalsPar2,method="BY");
pvalsHyb2_Par.adj <- p.adjust(pvalsHyb2_Par,method="BH")

# make output file
data.output2 <- cbind(data2[,1:5],(data2[,4]+data2[,5]),tempP2,pvalsPar2.adj,tempP2w, data2[,6:7], (data2[,6]+data2[,7]),tempH2,pvalsHyb2.adj,tempH2w,tempH2_P,pvalsHyb2_Par.adj);
colnames(data.output2) <- c("chrom", "start", "end", "ZHR", "Z30", "ZHR_Z30_total", "Parental.BETEst","Parental.LB","Parental.UB","Parental.BETP","Parental.BETQ","Parental.LB2","Parental.UB2","HYB_ZHR","HYB_Z30","HYB_TOTAL","HYB.BETEst","HYB.LB","HHYB.UB","HYB.BETP","HYB.BETQ","HYB.LBw","HYB.UBw","HYB_P.fetEst","HYB_P.LB","HYB_P.UB","HYB_P.P","HYB_P.Q")
write.table(data.output2,file = "ZHR_Z30_cis_trans_stats" ,sep="\t",quote=FALSE,row.names=FALSE)

################# cis and trans analysis/classification #####################

## make hyb ratios
hyb2_ratio<-(data2$HYB_ZHR/data2$HYB_Z30)
Par2_ratio<-(data2$ZHR/data2$Z30)

## log transform
log_hyb2_ratio <-log2(hyb2_ratio)
log_Par2_ratio <-log2(Par2_ratio)

## Make par/hyb ratio to properly determine cis+trans and cisXtrans
logPar2_loghyb2_rat <- log_Par2_ratio/log_hyb2_ratio

## Ratios to do analysis
ratios2 <- NULL
ratios2<-cbind(data.output2[1:5], data.output2[14:15], log_Par2_ratio, log_hyb2_ratio, data.output2$HYB.BETQ, data.output2$Parental.BETQ, data.output2$HYB_P.Q, logPar2_loghyb2_rat, pvalsHyb2, pvalsPar2, pvalsHyb2_Par)

# choose significant threshold
sig = 0.01

# Perform tests for determining classifications for each gene
a4<-subset(ratios2, ratios2[,10] > sig & ratios2[,11] > sig & ratios2[,12] > sig)
b4<-subset(ratios2, ratios2[,10] < sig & ratios2[,11] < sig & ratios2[,12] > sig)
c4<-subset(ratios2, ratios2[,10] > sig & ratios2[,11] < sig & ratios2[,12] < sig)
d4<-subset(ratios2, ratios2[,10] > sig & ratios2[,11] > sig & ratios2[,12] < sig)
e4<-subset(ratios2, ratios2[,10] < sig & ratios2[,11] > sig & ratios2[,12] > sig)
f4<-subset(ratios2, ratios2[,10] > sig & ratios2[,11] < sig & ratios2[,12] > sig)
g4<-subset(ratios2, ratios2[,10] < sig & ratios2[,11] > sig & ratios2[,12] < sig)
#cis+trans
h4<-subset(ratios2, ratios2[,10] < sig & ratios2[,11] < sig & ratios2[,12] < sig & ratios2[,8] > 0 &  ratios2[,9] > 0 & ratios2[,13] > 1)
i4<-subset(ratios2, ratios2[,10] < sig & ratios2[,11] < sig & ratios2[,12] < sig & ratios2[,8] < 0 &  ratios2[,9] < 0 & ratios2[,13] > 1)
#cisXtrans
j4<-subset(ratios2, ratios2[,10] < sig & ratios2[,11] < sig & ratios2[,12] < sig & ratios2[,8] > 0 &  ratios2[,9] < 0)
k4<-subset(ratios2, ratios2[,10] < sig & ratios2[,11] < sig & ratios2[,12] < sig & ratios2[,8] < 0 &  ratios2[,9] > 0)
l4<-subset(ratios2, ratios2[,10] < sig & ratios2[,11] < sig & ratios2[,12] < sig & ratios2[,8] > 0 &  ratios2[,9] > 0 & ratios2[,13] < 1)
m4<-subset(ratios2, ratios2[,10] < sig & ratios2[,11] < sig & ratios2[,12] < sig & ratios2[,8] < 0 &  ratios2[,9] < 0 & ratios2[,13] < 1)

var2<-c(nrow(a4),nrow(b4),nrow(c4),nrow(g4),(nrow(h4)+nrow(i4)),(nrow(j4)+nrow(k4)+nrow(l4)+nrow(m4)),(nrow(d4)+nrow(e4)+nrow(f4)),nrow(ratios2))
data.output5 <- matrix(data=var2,nrow=1, ncol=8, byrow=FALSE)
colnames(data.output5)<-c("Conserved","All cis","All trans","Compensatory","cis+trans","cisXtrans","Ambiguous","total")

## write tables and organize with classifcations
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_conserved <- as.data.frame(a4)
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_conserved$classification <- c("conserved")
colnames(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_conserved) <- c("chrom", "start", "end", "ZHR", "Z30", "HYB_ZHR", "HYB_Z30", "log_par_ratio", "log_hyb_ratio", "hyb_BETQ", "parental_BETQ", "hyb_p.q", "log_par_log_hyb_ratio", "hyb_p_vals", "par_p_vals", "hyb_par_p_vals", "classification")
write.table(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_conserved, file = "ZHR_Z30_reps_merged_downsamp_great20_0.05_conserved.txt")

ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis <- as.data.frame(b4)
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis$classification <- c("cis")
colnames(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis) <- c("chrom", "start", "end", "ZHR", "Z30", "HYB_ZHR", "HYB_Z30", "log_par_ratio", "log_hyb_ratio", "hyb_BETQ", "parental_BETQ", "hyb_p.q", "log_par_log_hyb_ratio", "hyb_p_vals", "par_p_vals", "hyb_par_p_vals", "classification")
write.table(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis, file = "ZHR_Z30_reps_merged_downsamp_great20_0.05_cis.txt")

ZHR_Z30_reps_merged_dm3_final_downsamp_great20_trans <- as.data.frame(c4)
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_trans$classification <- c("trans")
colnames(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_trans) <- c("chrom", "start", "end", "ZHR", "Z30", "HYB_ZHR", "HYB_Z30", "log_par_ratio", "log_hyb_ratio", "hyb_BETQ", "parental_BETQ", "hyb_p.q", "log_par_log_hyb_ratio", "hyb_p_vals", "par_p_vals", "hyb_par_p_vals", "classification")
write.table(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_trans, file = "ZHR_Z30_reps_merged_downsamp_great20_0.05_trans.txt")

ZHR_Z30_reps_merged_dm3_final_downsamp_great20_compens <- as.data.frame(g4)
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_compens$classification <- c("compensatory")
colnames(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_compens) <- c("chrom", "start", "end", "ZHR", "Z30", "HYB_ZHR", "HYB_Z30", "log_par_ratio", "log_hyb_ratio", "hyb_BETQ", "parental_BETQ", "hyb_p.q", "log_par_log_hyb_ratio", "hyb_p_vals", "par_p_vals", "hyb_par_p_vals", "classification")
write.table(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_compens, file = "ZHR_Z30_reps_merged_downsamp_great20_0.05_compensatory.txt")

ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis_p_trans <- as.data.frame(rbind(h4, i4))
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis_p_trans$classification <- c("cis_plus_trans")
colnames(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis_p_trans) <- c("chrom", "start", "end", "ZHR", "Z30", "HYB_ZHR", "HYB_Z30", "log_par_ratio", "log_hyb_ratio", "hyb_BETQ", "parental_BETQ", "hyb_p.q", "log_par_log_hyb_ratio", "hyb_p_vals", "par_p_vals", "hyb_par_p_vals", "classification")
write.table(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis_p_trans, file = "ZHR_Z30_reps_merged_downsamp_great20_0.05_cis_p_trans.txt")

ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis_x_trans <- as.data.frame(rbind(j4, k4, l4, m4))
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis_x_trans$classification <- c("cis_x_trans")
colnames(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis_x_trans) <- c("chrom", "start", "end", "ZHR", "Z30", "HYB_ZHR", "HYB_Z30", "log_par_ratio", "log_hyb_ratio", "hyb_BETQ", "parental_BETQ", "hyb_p.q", "log_par_log_hyb_ratio", "hyb_p_vals", "par_p_vals", "hyb_par_p_vals", "classification")
write.table(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis_x_trans, file = "ZHR_Z30_reps_merged_downsamp_great20_0.05_cis_x_trans.txt")

ZHR_Z30_reps_merged_dm3_final_downsamp_great20_ambig <- as.data.frame(rbind(d4, e4, f4))
ZHR_Z30_reps_merged_dm3_final_downsamp_great20_ambig$classification <- c("ambiguous")
colnames(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_ambig) <- c("chrom", "start", "end", "ZHR", "Z30", "HYB_ZHR", "HYB_Z30", "log_par_ratio", "log_hyb_ratio", "hyb_BETQ", "parental_BETQ", "hyb_p.q", "log_par_log_hyb_ratio", "hyb_p_vals", "par_p_vals", "hyb_par_p_vals", "classification")
write.table(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_ambig, file = "ZHR_Z30_reps_merged_downsamp_great20_0.05_ambig.txt")

###################### combine files to make master #####################

ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification <- rbind(as.data.frame(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_conserved), as.data.frame(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis), as.data.frame(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_trans), as.data.frame(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_compens), as.data.frame(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis_p_trans), as.data.frame(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_cis_x_trans), as.data.frame(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_ambig))
write.table(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification, file = "ZHR_Z30_reps_merged_downsamp_great20_0.05_all.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

####################  plot p val distribution ###########################

library(ggplot2)

ggplot(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification) +
  geom_density(aes(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification$hyb_p_vals))

ggplot(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification) +
  geom_density(aes(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification$par_p_vals))

ggplot(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification) +
  geom_density(aes(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification$hyb_par_p_vals))

######## plots ###########

ggplot(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification) +
  geom_point(aes(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification$log_par_ratio, ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification$log_hyb_ratio, color=classification, alpha=0.1)) +
  labs(y= "log(ZHR/Z30) HYB ALLELES", x = "log(ZHR/Z30) PARENTAL ALLELES") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept=c(-0,0), linetype="dotted") +
  geom_hline(yintercept=c(-0,0), linetype="dotted") +
  geom_abline(slope = 1, linetype="dotted")

## facet pot
ggplot(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification) +
  geom_point(aes(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification$log_par_ratio, ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification$log_hyb_ratio, alpha = 0.5, color=classification)) +
  labs(y= "log(ZHR/Z30) HYB ALLELES", x = "log(ZHR/Z30) PARENTAL ALLELES") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_vline(xintercept=c(-0,0), linetype="dotted") +
  geom_hline(yintercept=c(-0,0), linetype="dotted") +
  geom_abline(slope = 1, linetype="dotted") +
  facet_wrap( ~ classification, ncol=2)

######### plot effect sizes for all classes ##########
## parental
ggplot(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification, aes(x=classification, y=abs(log_par_ratio))) +
  geom_violin() +
  labs(y= "log2(ZHR/Z30) PARENTAL") +
  geom_boxplot(width=0.1)

## allelic
ggplot(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification, aes(x=classification, y=abs(log_hyb_ratio))) +
  geom_violin() +
  labs(y= "log2(ZHR/Z30) HYB ALLELES") +
  geom_boxplot(width=0.1)

### plot direction effect sizes
ggplot(ZHR_Z30_reps_merged_dm3_final_downsamp_great20_classification, aes(x=classification, y=log_par_ratio)) +
  geom_violin() +
  labs(y= "log2(ZHR/Z30) PARENTAL") +
  geom_boxplot(width=0.1)
