#!/usr/bin/env Rscript

library(pafr)

args = commandArgs(trailingOnly=TRUE)

ali<-read_paf(file.path(args[1]))
prefix<-file.path(args[2])

#alignment cov global
p1<-ggplot(ali, aes(alen, dv)) + 
    geom_point(alpha=0.6, colour="steelblue", size=2) + 
    scale_x_continuous("Alignment length (kb)", label =  function(x) x/ 1e3) +
    scale_y_continuous("Per base divergence") +
    theme_bw()

f1<-file.path(prefix,"assemblies.cov.pdf")
ggsave(f1)


#dot plot global
p2<-dotplot(ali,order_by="qstart") + theme_bw()
f2<-file.path(prefix,"assemblies.dot.pdf")
ggsave(f2)


#chrom_by_chrom
for (chrom in unique(ali$qname)) {

    to_keep <- list(
        c(chrom),
        c(agrep(chrom,unique(ali$tname),max = list(all = 8),value=TRUE))
    )


    p3<-dotplot(ali, label_seqs=TRUE, order_by="provided", ordering=to_keep)
    f3<-file.path(prefix, paste0(chrom ,".assemblies.dot.pdf"))
    ggsave(f3)

}

#maybe add others when we have sv calls