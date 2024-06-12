library(RSQLite)
library(stringr)
library(data.table)
library(TwoSampleMR)

####
#### take an example of one tissue "ctimp_Adipose_Subcutaneous.db"
####

########### find eQTLs on gene ESR1
filename = paste0("ctimp_Adipose_Subcutaneous.db")

sqlite.driver=dbDriver("SQLite")
db=dbConnect(sqlite.driver,dbname=filename)
dbListTables(db)
mytable=dbReadTable(db,"weights")
mytable[c('chr', 'pos','allele')] <- str_split_fixed(mytable$varID, '_', 3)
esr1 = mytable[mytable$chr == "chr6" & mytable$pos > 151640495 & mytable$pos < 152153274,]

# esr1 = mytable[mytable$gene == "ENSG00000091831.22",]
write.table(esr1,paste0("ctimp_Adipose_Subcutaneous.db.txt"),row.names = F,quote = F,col.names = T)

########### clump eQTLs
wt = as.data.frame(fread("ctimp_Adipose_Subcutaneous.db.txt",header=T))
colnames(wt)[c(2,7,8)] = c("SNP","chr_name","chrom_start")
wt$chr_name = gsub("chr","",wt$chr_name)

clp = clump_data(wt, clump_kb = 10000, clump_r2 = 0.2, clump_p1 = 1, clump_p2 = 1, pop = "EUR")
colnames(clp)[c(2,7,8)] = c("rsid","chr","pos")

write.table(clp,paste0("clump_ctimp_Adipose_Subcutaneous.db.txt"),row.names = F,quote = F,col.names = T)

############# make inference on causal effects

##disease
gdm = as.data.frame(fread("pheno_quantitative_40.csv",header=T))

##expression
traw = as.data.frame(fread("ESR1.traw",header = T)) ##indiviual-level genotype matrix on ESR1

eqtl = as.data.frame(fread(paste0("clump_ctimp_Adipose_Subcutaneous.db.txt"),header = T))

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
pv = out[2,4]
