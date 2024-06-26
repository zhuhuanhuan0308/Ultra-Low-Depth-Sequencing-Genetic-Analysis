library(data.table)
dat = as.data.frame(fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",header=T))
list = read.table("gene.list.txt",header=F)[,1]
sub = dat[dat$Description %in% list,]
sub2 = sub[,-c(1:2)]
mapply = apply(sub2,1,{function(x) (x-min(x))/(max(x)-min(x))})
# mapply = apply(sub2,1,{function(x) qqnorm(x,plot.it=F)$x})

slt = mapply[,-8]
rmeans = rowMeans(slt)
com = cbind.data.frame(colnames(sub2),rmeans)
out = com[order(com[,2],decreasing=T),]