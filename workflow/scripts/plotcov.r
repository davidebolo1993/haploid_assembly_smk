#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

df<-fread(file.path(args[1]))
pdf<-file.path(args[2])
df<-subset(df, V1 == "total")

p<-ggplot(df, aes(x=V2, y=V3,col=V5)) + 
    geom_line() + 
    labs(x="Depth", y=expression("Fraction of bases ">=" depth")) + 
    theme_bw() + 
    theme(legend.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal") + 
    scale_x_continuous(limits=c(0,200)) +
    facet_wrap(~V4, nrow=1)

ggsave(pdf, height=10, width=15)