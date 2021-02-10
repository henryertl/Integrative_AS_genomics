# reorder yeast data

# read in matrices
rep1 <- read.delim(file="/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/evx035f3p_Supp/Rep_1.txt", header = T) %>% as.data.frame()
rep2 <- read.delim(file="/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/evx035f3p_Supp/Rep_2.txt", header = T) %>% as.data.frame()
rep1_2 <- merge(rep1, rep2, by = 'Gene') %>% as.data.frame()

# CPM transform whole matrix
rep1_2_CPM <- matrix(ncol = ncol(rep1_2) - 1, nrow = nrow(rep1_2)) %>% as.data.frame()

for(i in 2:ncol(rep1_2)) {
  rep1_2_CPM[,i] <- (rep1_2[,i]/sum(rep1_2[,i]))*1000000
}

# reassign headers and gene namees
colnames(rep1_2_CPM) <- colnames(rep1_2)
rep1_2_CPM$Gene <- rep1_2$Gene

# rearrange to have replicates together and strain comparisons separate
C_R <- rep1_2_CPM[, c("Gene", "C1.x", "C1.y", "R.x", "R.y", "CxR.C.x", "CxR.C.y", "CxR.R.x", "CxR.R.y")]
C_P <- rep1_2_CPM[, c("Gene", "C2.x", "C2.y", "P.x", "P.y", "CxP.C.x", "CxP.C.y", "CxP.P.x", "CxP.P.y")]
C_M <- rep1_2_CPM[, c("Gene", "C3.x", "C3.y", "M.x", "M.y", "CxM.C.x", "CxM.C.y", "CxM.M.x", "CxM.M.y")]
C_B <- rep1_2_CPM[, c("Gene", "C4.x", "C4.y", "B.x", "B.y", "CxB.C.x", "CxB.C.y", "CxB.B.x", "CxB.B.y")]

# write files to save
write.table(rep1_2_CPM, file = "/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/four_yeast_AS_RNA_CPM_ALL.txt", row.names = F, quote = F, sep = "\t")
write.table(C_R, file = "/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/four_yeast_AS_RNA_CPM_C_R.txt", row.names = F, quote = F, sep = "\t")
write.table(C_P, file = "/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/four_yeast_AS_RNA_CPM_C_M.txt", row.names = F, quote = F, sep = "\t")
write.table(C_M, file = "/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/four_yeast_AS_RNA_CPM_C_P.txt", row.names = F, quote = F, sep = "\t")
write.table(C_B, file = "/Users/henryertl/Documents/Wittkopp_lab/Four_spp_yeast/four_yeast_AS_RNA_CPM_C_B.txt", row.names = F, quote = F, sep = "\t")
