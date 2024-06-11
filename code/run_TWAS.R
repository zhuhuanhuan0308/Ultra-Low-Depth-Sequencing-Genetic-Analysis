args = commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])

library(RSQLite)
library(stringr)

files = read.table("/hwfssz5/ST_HEALTH/P18Z10200N0124/zhuhuanhuan1/NIPT/Wuhan/MR/TWAS/ctimp.idx.txt",header = F)[,1]
filename = paste0("/zfssz3/pub/database/zenodo.org/records/5709385/files/",files[i])

sqlite.driver=dbDriver("SQLite")
db=dbConnect(sqlite.driver,dbname=filename)
dbListTables(db)
mytable=dbReadTable(db,"weights")
mytable[c('chr', 'pos','allele')] <- str_split_fixed(mytable$varID, '_', 3)
esr1 = mytable[mytable$chr == "chr6" & mytable$pos > 151640495 & mytable$pos < 152153274,]

# esr1 = mytable[mytable$gene == "ENSG00000091831.22",]
write.table(esr1,paste0("/hwfssz5/ST_HEALTH/P18Z10200N0124/zhuhuanhuan1/NIPT/Wuhan/MR/TWAS/ESR1_CTIMP/",files[i],".txt"),row.names = F,quote = F,col.names = T)

###############################################################
library(data.table)
library(TwoSampleMR)

setwd("C:/Users/zhuhuanhuan1/Downloads/BGI/Projects/NIPT/community review/revise/GDM/MR/TWAS/ESR1_CTIMP")

files = list.files(".",full.names = F)
nf = length(files)

for(i in 1:nf){
  wt = as.data.frame(fread(files[i],header=T))
  colnames(wt)[c(2,7,8)] = c("SNP","chr_name","chrom_start")
  wt$chr_name = gsub("chr","",wt$chr_name)
  
  clp = clump_data(wt, clump_kb = 10000, clump_r2 = 0.2, clump_p1 = 1, clump_p2 = 1, pop = "EUR")
  colnames(clp)[c(2,7,8)] = c("rsid","chr","pos")
  
  write.table(clp,paste0("C:/Users/zhuhuanhuan1/Downloads/BGI/Projects/NIPT/community review/revise/GDM/MR/TWAS/clump_ESR1_CTIMP/",files[i]),row.names = F,quote = F,col.names = T)
}


#################################################################

library(data.table)

#disease
gdm = as.data.frame(fread("/hwfssz5/ST_HEALTH/P18Z10200N0124/lilinxuan/workshop/20220519_wuhan_phenotype_final/phenotype_genaral/pheno_quantitative_40.csv",header=T))

##  /hwfssz5/ST_HEALTH/P18Z10200N0124/zhuhuanhuan1/NIPT/Wuhan/MR/pheno_GDM.csv

##expression
traw = as.data.frame(fread("/hwfssz5/ST_HEALTH/P18Z10200N0124/lilinxuan/workshop/20231121_extract_snps/output/ESR1.traw",header = T))

files = read.table("/hwfssz5/ST_HEALTH/P18Z10200N0124/zhuhuanhuan1/NIPT/Wuhan/MR/TWAS/ctimp.idx.txt",header = F)[,1]

nf = length(files)
pv = rep(NA,nf)

for(i in 1:nf){
  eqtl = as.data.frame(fread(paste0("/hwfssz5/ST_HEALTH/P18Z10200N0124/zhuhuanhuan1/NIPT/Wuhan/MR/TWAS/ESR1_CTIMP/",files[i],".txt"),header = T))
  
  snp = intersect(eqtl$pos, traw$POS)
  
  subeqtl = eqtl[match(snp,eqtl$pos),]
  subtraw = traw[match(snp,traw$POS),]
  print(identical(subeqtl$ref_allele,subtraw$COUNTED))
  
  subeqtl[which(subeqtl$ref_allele!=subtraw$COUNTED),"weight"] = -subeqtl[which(subeqtl$ref_allele!=subtraw$COUNTED),"weight"]
  
  prs = t(subtraw[,-c(1:6)]) %*% (subeqtl$weight)
  
  ## disease
  smp = intersect(row.names(prs),paste0("0_",gdm$IID))
  
  subprs = prs[match(smp,row.names(prs)),]
  subgdm = gdm[match(smp,paste0("0_",gdm$IID)),]
  
  out = summary(glm(subgdm$OGLU~subprs,family = "gaussian"))$coefficients
  pv[i] = out[2,4]
  
  write.table(out,paste0("/hwfssz5/ST_HEALTH/P18Z10200N0124/zhuhuanhuan1/NIPT/Wuhan/MR/TWAS/CTIMP/OGLU_",files[i],".txt"),row.names = T,quote = F,col.names = T)
  
}

pvs = cbind.data.frame(files,pv)
pvs[order(pvs$pv),]
# write.table(pvs,"/hwfssz5/ST_HEALTH/P18Z10200N0124/zhuhuanhuan1/NIPT/Wuhan/MR/TWAS/CTIMP.txt",row.names=F,quote=F,col.names=T)