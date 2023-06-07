#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

df<-fread(file.path(args[1]))
pdf<-file.path(args[2])

p<-ggplot(data=df, aes(x=V3, y=V4)) + 
  stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  scale_fill_distiller(palette=4, direction=1) +
  scale_x_continuous("read_length (bp)",expand=c(0,0)) + 
  scale_y_continuous("mean_quality (phred)",expand=c(0,0)) + 
  theme_bw()+
  guides(alpha="none",col=guide_legend(title="Density"))+
  theme(legend.position='right', legend.background=element_blank(),legend.direction="vertical")+
  facet_wrap(~V8, scales = "fixed", nrow=1)+
  guides(fill=guide_legend(title="Density"))

ggsave(pdf, height=10, width=15)
