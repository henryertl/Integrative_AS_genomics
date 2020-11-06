###############################################################################
################# Coolon Method for ASE cis - trans analysis ##################
###############################################################################

######################### hypergeometric sampling #############################

# read in datasets
ZHR_Z30_all <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Data_tables/Analysis_regions_files/ZHR_Z30_ATAC_merged_counts_final_peakset.bed", header = F)
ZHR_TSIM_all <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Data_tables/Analysis_regions_files/ZHR_TSIM_ATAC_merged_counts_final_peakset.bed", header = F)
ZHR_Z30_TSIM_all <- cbind(ZHR_Z30_all, ZHR_TSIM_all)
ZHR_Z30_TSIM_all <- ZHR_Z30_TSIM_all[c(1:15,19:ncol(ZHR_Z30_TSIM_all))]
colnames(ZHR_Z30_TSIM_all) <- c("chrom", "start", "end", "ZHR_1", "ZHR_2", "ZHR_3", "Z30_1", "Z30_2", "Z30_3", "HYB_1_ZHR", "HYB_2_ZHR", "HYB_3_ZHR", "HYB_1_Z30", "HYB_2_Z30", "HYB_3_Z30", "ZHR_1_tsim", "ZHR_2_tsim", "ZHR_3_tsim", "TSIM_1", "TSIM_2", "TSIM_3", "HYB_1_ZHR_tsim", "HYB_2_ZHR_tsim", "HYB_3_ZHR_tsim", "HYB_1_TSIM", "HYB_2_TSIM", "HYB_3_TSIM")



# get total read counts
sum(full_dataset[,4])
# 367122
sum(full_dataset[,5])
# 776279
sum(full_dataset[,6])
# 625190
sum(full_dataset[,7])
min <- 227904
sum(full_dataset[,8])
# 545868
sum(full_dataset[,9])
# 521377
sum(full_dataset[,10])
# 595665
sum(full_dataset[,11])
# 553832

ZHR_Z30_reps_merged_dm3_final <- NULL
ZHR_Z30_reps_merged_dm3_final <- full_dataset %>% as.data.frame()

# randomly down sample each to 227904 reads and add to new file
ZHR_Z30_reps_merged_dm3_final$number <- (1:nrow(ZHR_Z30_reps_merged_dm3_final))

ZHR_Z30_reps_merged_dm3_final_downsamp <- NULL %>% as.data.frame()
ZHR_Z30_reps_merged_dm3_final_downsamp <- ZHR_Z30_reps_merged_dm3_final[c(1:3)]
ZHR_Z30_reps_merged_dm3_final_downsamp$ID <- cbind(1:nrow(ZHR_Z30_reps_merged_dm3_final)) %>% as.data.frame()

# initiate variables



######### ZHR  1 ############
# expand read counts and sample by min number
a1<-c()
for (i in ZHR_Z30_reps_merged_dm3_final[,12]) {
  temp<-rep.int(i,ZHR_Z30_reps_merged_dm3_final[i,4])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "ZHR_1_downsamp")
d_merge <- merge(d_table, ZHR_Z30_reps_merged_dm3_final, by.x = "ID", by.y = "number")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)

######### ZHR  3 ############
a1<-c()
for (i in ZHR_Z30_reps_merged_dm3_final[,12]) {
  temp<-rep.int(i,ZHR_Z30_reps_merged_dm3_final[i,5])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "ZHR_3_downsamp")
d_merge2 <- merge(d_table, d_merge, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge)

######### Z30  1 ############
a1<-c()
for (i in ZHR_Z30_reps_merged_dm3_final[,12]) {
  temp<-rep.int(i,ZHR_Z30_reps_merged_dm3_final[i,6])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "Z30_1_downsamp")
d_merge3 <- merge(d_table, d_merge2, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge2)

######### Z30  2 ############
a1<-c()
for (i in ZHR_Z30_reps_merged_dm3_final[,12]) {
  temp<-rep.int(i,ZHR_Z30_reps_merged_dm3_final[i,7])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "Z30_2_downsamp")
d_merge4 <- merge(d_table, d_merge3, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge3)

######### HYB_1_ZHR ############
a1<-c()
for (i in ZHR_Z30_reps_merged_dm3_final[,12]) {
  temp<-rep.int(i,ZHR_Z30_reps_merged_dm3_final[i,8])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_1_ZHR_downsamp")
d_merge5 <- merge(d_table, d_merge4, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge4)

######### HYB_1_Z30 ############
a1<-c()
for (i in ZHR_Z30_reps_merged_dm3_final[,12]) {
  temp<-rep.int(i,ZHR_Z30_reps_merged_dm3_final[i,9])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_1_Z30_downsamp")
d_merge6 <- merge(d_table, d_merge5, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge5)

######### HYB_2_ZHR ############
a1<-c()
for (i in ZHR_Z30_reps_merged_dm3_final[,12]) {
  temp<-rep.int(i,ZHR_Z30_reps_merged_dm3_final[i,10])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_2_ZHR_downsamp")
d_merge7 <- merge(d_table, d_merge6, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge6)

######### HYB_2_Z30 ############
a1<-c()
for (i in ZHR_Z30_reps_merged_dm3_final[,12]) {
  temp<-rep.int(i,ZHR_Z30_reps_merged_dm3_final[i,11])
  a1<-c(a1,temp)
}

d<-sample(a1,min,replace = FALSE)

# collapse, count, and append to new file
d_table <- table(d)
d_table <- as.data.frame(d_table)
colnames(d_table) <- c("ID", "HYB_2_Z30_downsamp")
d_merge8 <- merge(d_table, d_merge7, by.x = "ID", by.y = "ID")

rm(d)
rm(d_table)
rm(dt)
rm(a1)
rm(c)
rm(i)
rm(temp)
rm(d_merge7)

d_merge8 <- as.data.frame(d_merge8)
write.table(d_merge8, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/ZHR_Z30_reps_merged_dm3_final_downsamp_all_ALL.txt", sep = "\t", row.names = F, quote = F)
down_samp_only <- NULL
down_samp_only <- d_merge8[c(10:12,9,8,7,6,5,3,4,2)]
write.table(down_samp_only, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/ZHR_Z30_reps_merged_dm3_final_downsamp_ONLY.txt", sep = "\t", row.names = F, quote = F)



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
