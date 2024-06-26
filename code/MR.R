library(TwoSampleMR)
library(data.table)
library(R.utils)

ss = read.table("file_list.txt",header = T)[,1]
cc = as.data.frame(fread("case-control-2.txt",header = T))
idx = read.table("idx.txt",header = F)[,1]

####
which(pvals < 0.05/(51*46),arr.ind = T)

cl = cc[,2]
which(cl == "Asthma")
which(cl == 'Graves_disease')
which(cl == 'Prostate_cancer')
which(cl == 'Urolithiasis')
####
n_wh = length(ss)
n_bbj = nrow(cc)

pvals = read.table('MR.pvals.final.txt',header = T)

ord <- fread('catagory_order.txt')
pvals <- pvals[as.numeric(na.omit(match(ord$Phenoname, rownames(pvals)))),]
ord <- ord[as.numeric(na.omit(match(rownames(pvals),ord$Phenoname))),]

rownames(pvals) <- ord$Abbreviation

colnames(pvals) <- c("Rheumatoid arthritis in Asian",
                     "Open-angle glaucoma",
                     "Type 2 Diabetes",
                     "Smoking initiation",
                     "Smoking cessation",
                     "Arrhythmia",
                     "Asthma",
                     "Atopic dermatitis",
                     "Biliary tract cancer",
                     "Cataract",
                     "Cerebral aneurysm",
                     "Cervical cancer",
                     "Chronic hepatitis B",
                     "Chronic hepatitis C",
                     "Chronic obstructive pulmonary disease",
                     "Cirrhosis",
                     "Congestive heart failure",
                     "Drug eruption",
                     "Endometrial cancer",
                     "Endometriosis",
                     "Epilepsy",
                     "Esophageal cancer",
                     "Gastric cancer",
                     "Glaucoma",
                     "Graves disease",
                     "Hematological malignancy",
                     "Interstitial lung disease",
                     "Ischemic stroke",
                     "Keloid",
                     "Lung cancer",
                     "Nephrotic syndrome",
                     "Osteoporosis",
                     "Ovarian cancer",
                     "Pancreatic cancer",
                     "Periodontal disease",
                     "Peripheral artery disease",
                     "Pollinosis",
                     "Prostate cancer",
                     "Pulmonary tuberculosis",
                     "Rheumatoid arthritis",
                     "Type 2 diabetes",
                     "Urolithiasis",
                     "Uterine fibroids",
                     "hepatocellular carcinoma",
                     "Breast cancer",
                     "Coronary artery disease")
  
  dd <- matrix(data = c(35,5,4,17,49,10,29,41,42,42), byrow = F, ncol = 2)
for (x in 1:5) {
  i = dd[x,1]
  j = dd[x,2]

  wh = as.data.frame(fread(paste0("WH/",ss[i],".ss"),header = T))
  wh$pvalue = as.numeric(wh$pvalue)
  wh_sig = wh[wh$pvalue < 5e-08 & wh$SNP != ".",]
  exposure_dat = format_data(
    wh_sig,
    type = "exposure",
    phenotype_col = "trait",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pvalue",
    samplesize_col = "samplesize"
  )
  exposure_dat$chr.exposure = gsub("chr","",exposure_dat$chr.exposure)
  
  exp_dat = clump_data(
    exposure_dat,
    clump_kb = 10000,
    clump_r2 = 0.1,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EAS"
  )
  

  trait = cc[j,2]
  ncase = cc[j,3]
  ncontrol = cc[j,4]
  
  bbj = as.data.frame(fread(paste0("BBJ/","phenocode-",trait,".tsv.gz"),header = T))
  bbj_slt = bbj[bbj$rsids %in% exp_dat$SNP,]
  bbj_slt$trait = trait
  bbj_slt$ncase = ncase
  bbj_slt$ncontrol = ncontrol
  
  out_dat = format_data(
    bbj_slt,
    type = "outcome",
    phenotype_col = "trait",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "maf",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    ncase_col = "ncase",
    ncontrol_col = "ncontrol",
  )
  
  dat <- harmonise_data(exp_dat, out_dat)
  res <- mr(dat)

  mr_pleiotropy_test(dat)
  mr_heterogeneity(dat)
  
  pdf(paste0(ss[i],'_VS_',trait,'.pdf'),height = 4,width = 4)
  single <- mr_leaveoneout(dat)
  mr_leaveoneout_plot(single)
  dev.off()
  x = x+1
}  
  
library(MRPRESSO)
library(MendelianRandomization)

mr_scatter_plot(res, dat)
mr_plot(MRAllObject_all)


dd <- matrix(data = c(35,5,17,49,10,29,42,46), byrow = F, ncol = 2)
a <- matrix(NA, ncol = 11, nrow = 4)
for (x in c(1:4)) {
  i = dd[x,1]
  #35,5,4,17,49
  j = dd[x,2]
  #10,29,41,42,42
  wh = as.data.frame(fread(paste0("WH/",ss[i],".ss"),header = T))
  wh$pvalue = as.numeric(wh$pvalue)
  wh_sig = wh[wh$pvalue < 5e-08 & wh$SNP != ".",]
  exposure_dat = format_data(
    wh_sig,
    type = "exposure",
    phenotype_col = "trait",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pvalue",
    samplesize_col = "samplesize"
  )
  exposure_dat$chr.exposure = gsub("chr","",exposure_dat$chr.exposure)
  
  exp_dat = clump_data(
    exposure_dat,
    clump_kb = 10000,
    clump_r2 = 0.1,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EAS"
  )
  trait = cc[j,2]
  ncase = cc[j,3]
  ncontrol = cc[j,4]
  
  bbj = as.data.frame(fread(paste0("BBJ/","phenocode-",trait,".tsv.gz"),header = T))
  bbj_slt = bbj[bbj$rsids %in% exp_dat$SNP,]
  bbj_slt$trait = trait
  bbj_slt$ncase = ncase
  bbj_slt$ncontrol = ncontrol
  
  out_dat = format_data(
    bbj_slt,
    type = "outcome",
    phenotype_col = "trait",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "maf",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    ncase_col = "ncase",
    ncontrol_col = "ncontrol",
  )
  
  dat <- harmonise_data(exp_dat, out_dat)
  res <- mr(dat)
  pleiotropy <- mr_pleiotropy_test(dat)
  ivw <- mr_ivw(dat$beta.exposure,dat$beta.outcome,dat$se.exposure,dat$se.outcome)
  heterogeneity <- mr_heterogeneity(dat)
  egger <- mr_egger_regression(dat$beta.exposure,dat$beta.outcome,dat$se.exposure,dat$se.outcome)
  
  a[x,1] <- ss[i]
  a[x,2] <- trait
  a[x,3] <- ivw$b
  a[x,4] <- ivw$se
  a[x,5] <- ivw$pval
  a[x,6] <- pleiotropy$egger_intercept
  a[x,7] <- pleiotropy$se
  a[x,8] <- pleiotropy$pval
  a[x,9] <- heterogeneity$Q[heterogeneity$method == 'Inverse variance weighted']
  a[x,10] <- heterogeneity$Q_pval[heterogeneity$method == 'Inverse variance weighted']
  a[x,11] <- egger$nsnp
}

colnames(a) <- c('Exposure','Outcome', 'Beta_IVW', 'SE_IVW', 'Pvalue_IVW', 'Intercept _pleiotropy', 'SE_pleiotropy', 'Pvalue_pleiotropy', 'Q', 'Q_Pvalue', 'nsnp')
write.table(a, file = 'result.txt', quote = F, col.names = T, row.names = F,sep = '\t')