

Grh_ATAC_all <- read.delim("/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/Grh_ATAC_ZHR_Z30.bed", header = F)

Grh_ATAC_all$V10 <- paste(Grh_ATAC_all$V1, Grh_ATAC_all$V2, Grh_ATAC_all$V3, sep="_")

ATAC_Grh <- subset(Grh_ATAC_all, Grh_ATAC_all$V10 != ".") %>% as.data.frame()
ATAC_no_Grh <- subset(Grh_ATAC_all, Grh_ATAC_all$V10 == ".") %>% as.data.frame()
Grh_ATAC_all_final <- Grh_ATAC_all[c(10,4,8)] %>% as.data.frame()


colnames(Grh_ATAC_all_final) <- c("Paste_locus", "ZHR_ATAC", "Grh_bind")

Grh_ATAC_all_final$Grh_bind <- Grh_ATAC_all_final$Grh_bind %>% as.numeric()

Grh_ATAC_all_counts <- ddply(Grh_ATAC_all_final,.(Paste_locus),summarize,ZHR_ATAC=sum(ZHR_ATAC), Grh_bind=sum(Grh_bind))

ggplot(Grh_ATAC_all_counts, aes(x=ZHR_ATAC, y=Grh_bind)) +
geom_point() +
xlim(0,1000) +
ylim(0,1000)

ggplot(Grh_ATAC_all_counts,aes(x=ZHR_ATAC)) +
geom_density() +
xlim(0,1000)

ggplot(Grh_ATAC_all_counts,aes(x=Grh_bind)) +
geom_histogram() +
xlim(0,1000)

ggplot(Grh_ATAC_all_final,aes(x=Grh_bind)) +
geom_histogram() +
xlim(0,1000)
