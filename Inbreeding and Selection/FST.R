library(emdbook)
library(ggplot2)
library(dplyr)
library(qqman)

#FST Plot
data<-read.csv("~/Desktop/FST/dogs_australia_fst.windowed.weir.fst", sep = "", header = T)

#Mean FST
data_100 <- data %>% filter(N_VARIANTS > 50)
top_0.1 <- quantile(data_100$MEAN_FST, 0.999)
mean(data_100$MEAN_FST)
manhattan(data_100, chr="CHROM", bp="BIN_START", snp = "BIN_START", p="MEAN_FST", logp=F, genomewideline=top_0.1, col = c("black", "grey"))

final<-data_100 %>% filter(MEAN_FST > top_0.1)

MEAN_FST_0.1 <- data_100 %>%
  filter(MEAN_FST > top_0.1)
write.csv(MEAN_FST_0.1, "~/Desktop/MEAN_FST_0.1.csv")

#Weighted FST
top_0.1 <- quantile(data_100$WEIGHTED_FST, 0.999)
mean(data_100$WEIGHTED_FST)
manhattan(data_100, chr="CHROM", bp="BIN_START", snp = "BIN_START", p="WEIGHTED_FST", logp=F, genomewideline=top_0.1, col = c("black", "grey"))

WEIGHTED_FST_0.1 <- data_100 %>%
  filter(WEIGHTED_FST > top_0.1)
write.csv(WEIGHTED_FST_0.1, "~/Desktop/WEIGHTED_FST_0.1.csv")

