# get ratios of reinforcing - opposing cis-trans and plot
library(dplyr)
library(ggplot2)

theme_main <- function() {
  theme_bw() +
  theme(
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 20),
  strip.text = element_text(size = 20),
  legend.text= element_text(size = 15),
  legend.title = element_text(size = 20),
  plot.title = element_text(size = 25, face = "bold")
)
}

## RNA_seq

df <- read.delim("~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/ZHR_Z30_TSIM_Full_results_output.txt", header = T)
x <- nrow(df[df$Direction.x == "Opposing",])
y <- nrow(df[df$Direction.x == "Reinforcing",])
intra_opposing <- x/(x+y)
intra_reinforcing <- y/(x+y)

a <- nrow(df[df$Direction.y == "Opposing",])
b <- nrow(df[df$Direction.y == "Reinforcing",])
inter_opposing <- a/(a+b)
inter_reinforcing <- b/(a+b)

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("Mel-Mel", 2), rep("Mel-Sim", 2))
twobytwo[,2] <- c(intra_opposing,intra_reinforcing,inter_opposing,inter_reinforcing)
twobytwo[,3] <- c("Cis-trans opposing","Cis-trans reinforcing","Cis-trans opposing","Cis-trans reinforcing")
colnames(twobytwo) <- c("Comparison", "Proportion", "Cis_trans_Directionality")

A <- ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Cis_trans_Directionality)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  ggtitle("Gene expression: Cis-trans opposing vs reinforcing within & between species")

ggsave(A, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/RNA_seq/Figures/cis_trans_opposing_reinforcing_RNA.pdf", width = 15, height = 15)

## ATAC_seq
df <- read.delim("~/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Bayes_test_outputs/ZHR_Z30_TSIM_universal_peakset_outputs/ZHR_Z30_TSIM_Full_results_output.txt", header = T)
x <- nrow(df[df$Direction.x == "Opposing",])
y <- nrow(df[df$Direction.x == "Reinforcing",])
intra_opposing <- x/(x+y)
intra_reinforcing <- y/(x+y)

a <- nrow(df[df$Direction.y == "Opposing",])
b <- nrow(df[df$Direction.y == "Reinforcing",])
inter_opposing <- a/(a+b)
inter_reinforcing <- b/(a+b)

twobytwo <- matrix(ncol = 3, nrow = 4) %>% as.data.frame()
twobytwo[,1] <- c(rep("Mel-Mel", 2), rep("Mel-Sim", 2))
twobytwo[,2] <- c(intra_opposing,intra_reinforcing,inter_opposing,inter_reinforcing)
twobytwo[,3] <- c("Cis-trans opposing","Cis-trans reinforcing","Cis-trans opposing","Cis-trans reinforcing")
colnames(twobytwo) <- c("Comparison", "Proportion", "Cis_trans_Directionality")

B <- ggplot(data = twobytwo, aes(x=Comparison, y=Proportion, fill=Cis_trans_Directionality)) +
  geom_bar(position="stack", stat="identity") +
  xlab("") +
  theme_main() +
  ggtitle("Chromatin Accessibility: Cis-trans opposing vs reinforcing within & between species")

ggsave(B, file = "/Users/henryertl/Documents/Wittkopp_lab/AS_ATAC_RNA_2020_10_1/ATAC_seq/Figures/cis_trans_opposing_reinforcing_ATAC.pdf", width = 15, height = 15)
