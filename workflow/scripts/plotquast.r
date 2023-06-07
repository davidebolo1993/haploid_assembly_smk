#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

df<-fread(file.path(args[1]))
pdf<-file.path(args[2])
vals<-c("# contigs (>= 0 bp)",  "Total length (>= 0 bp)", "NG50", "LG50", "Genome fraction (%)", "K-mer-based compl. (%)", "Complete BUSCO (%)", "Partial BUSCO (%)", "# misassemblies", "Misassembled contigs length")
df<-subset(df, V1 %in% vals)
df$V1<-factor(df$V1, levels=vals)
df$V2<-as.numeric(df$V2)


p<-ggplot(df, aes(x=V4,y=V2,fill=V3,width=.3)) +
  geom_bar(stat="identity", position="dodge")+
  facet_wrap(~V1, scales='free_y', nrow=1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.title=element_blank()) +
  labs(x="",y="")

ggsave(pdf, height=10, width=25)