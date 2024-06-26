# Ultra-Low-Depth Sequencing Genetic Analysis

## Overview
This repository serves as a dedicated space for housing the codebase used in our publication for genetic analysis of ultra-low depth sequencing data. It is intended to provide researchers with access to the methodologies and algorithms employed in our study, facilitating further research and analysis in the field of genetic.

## License
The code within this repository is licensed under the [MIT License](./LICENSE). Please refer to the license file for more information on the terms and conditions of using and contributing to this project.

For more detailed information on the methodology and results, please refer to our accompanying publication:
1. https://doi.org/10.1101/2023.11.23.23298979
2. https://doi.org/10.1101/2023.11.23.23298977

## 1. Sequecing data preprocessing and genotype imputaion
### 1.1. Alignment
code：workflow_bwa.sh
```bash
### static path for tools and reference on local disk.
### Tools
gatk=./software/gatk-4.0.4.0/gatk
java=./software/jdk1.8.0_131/bin/java
picard=./software/picard-2.10.10/picard.jar
bwa=./software/bwa/0.7.16/bwa
samtools=./software/samtools
fastp=./software/fastp-0.23.2/fastp

### Reference Genome
hg38=./toolset/reference/hg38_test/Homo_sapiens_assembly38.fasta

### GATK bundle
gatk_bundle_dir=./database/ftp.broadinstitute.org/gsapubftp-anonymous/bundle/hg38
known_indels=./pub/database/ftp.broadinstitute.org/gsapubftp-anonymous/bundle/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz

### Input

sample_id=$1
lane_id=$2
fq=$3
out_path=$4

### Output
final_outdir=$out_path/final_data/$sample_id
outdir=$out_path/temp_data/$sample_id

if [ ! -d $final_outdir ]
then mkdir -p $final_outdir
fi

if [ ! -d $outdir ]
then mkdir -p $outdir
fi

######################################################################################
################################### Pipeline #########################################
######################################################################################
echo "We're doing the job in $sample_id"
echo "We are calculating $fq"
echo "We are doing $lane_id"
echo "We'll save it in $final_outdir"

### step 0: fastp for QC

time $fastp -i $fq -o $outdir/${lane_id}.clean.fq.gz --qualified_quality_phred=5 --unqualified_percent_limit=50 --n_base_limit=10 \
--adapter_sequence="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA" --adapter_sequence_r2="AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG" \
--disable_trim_poly_g --thread=16 -j $outdir/${lane_id}.json -h $outdir/${lane_id}.html -R $lane_id

### step 1: bwa

echo ""
time $bwa aln -e 10 -t 4 -i 5 -q 0 $hg38 $outdir/${lane_id}.clean.fq.gz > $outdir/${lane_id}.sai && \
    time $bwa samse -r "@RG\tID:${sample_id}\tPL:COMPLETE\tSM:${lane_id}" $hg38 $outdir/${lane_id}.sai $outdir/${lane_id}.clean.fq.gz | $samtools view -h -Sb - > $outdir/${lane_id}.bam && echo "** bwa done **" && \
    time $samtools sort -@ 8 -O bam -o $outdir/${lane_id}.sorted.bam $outdir/${lane_id}.bam && echo "** bam sorted done **" && \
    time $samtools rmdup $outdir/${lane_id}.sorted.bam $outdir/${lane_id}.sorted.rmdup.bam && echo "** rmdup done **" && \
    time $samtools index $outdir/${lane_id}.sorted.rmdup.bam && echo "** index done **" && touch ${outdir}/bwa_sort_rmdup.finish

if [ ! -f ${outdir}/bwa_sort_rmdup.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bwa_sort_rmdup not done **" && exit
fi

### step 2: gatk

echo ""
time $gatk BaseRecalibrator \
    -R $hg38 \
    -I $outdir/${lane_id}.sorted.rmdup.bam \
    --known-sites $gatk_bundle_dir/dbsnp_146.hg38.vcf.gz \
    --known-sites $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites $known_indels \
    -O $outdir/${lane_id}.recal_data.table && echo "** BaseRecalibrator done " && touch ${outdir}/baseRecalibrator.finish

if [ ! -f ${outdir}/baseRecalibrator.finish ]
then echo "** [WORKFLOW_ERROR_INFO] baseRecalibrator not done **" && exit
fi

time $gatk ApplyBQSR \
    -R $hg38 \
    --bqsr-recal-file $outdir/${lane_id}.recal_data.table \
    -I $outdir/${lane_id}.sorted.rmdup.bam \
    -O $outdir/${lane_id}.sorted.rmdup.BQSR.bam && echo "** PrintReads done **" && touch ${outdir}/PrintReads.finish

if [ ! -f ${outdir}/PrintReads.finish ]
then echo "** [WORKFLOW_ERROR_INFO] PrintReads not done **" && exit
fi

### step 3: bam index

time $samtools index $outdir/${lane_id}.sorted.rmdup.BQSR.bam && echo "** bam index done **" && touch ${outdir}/bam_index.finish

if [ ! -f ${outdir}/bam_index.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bam index not done **" && exit
fi

### step 4: bam stats
time $samtools stats $outdir/${lane_id}.sorted.rmdup.BQSR.bam > $outdir/${lane_id}.sorted.rmdup.BQSR.bamstats && echo "** bamstats done **" && touch ${outdir}/bamstats.finish

if [ ! -f ${outdir}/bamstats.finish ]
then echo "** [WORKFLOW_ERROR_INFO] bamstats not done **" && exit
fi

### move to final_dir

mv -f $outdir/${lane_id}.sorted.rmdup.BQSR.bam* $final_outdir && echo "** move2final done **" && touch ${outdir}/move2final.finish
if [ ! -f ${outdir}/move2final.finish ]
then echo "** [WORKFLOW_ERROR_INFO] move2final not done **" && exit #v1.6
fi

### clear up
rm -vrf $outdir
```

usage, for example:
``` bash
bash workflow_bwa.sh CL100045254_L01_23 CL100045254_L01_23 CL100045254_L01_23.fq.gz output_path/
```

### 1.2. genotype imputation
code: STITCH.R
```R
#! R
args=commandArgs(T)
outputdir=args[1]
bamlist=args[2]
ref=args[3]
human_posfile=args[4]
human_K=as.numeric(args[5])
human_nGen=as.numeric(args[6])
nCores=as.numeric(args[7])
niterations=as.numeric(args[8])
regionStart=as.numeric(args[9])
regionEnd=as.numeric(args[10])
chr=as.character(args[11])
buffer=as.numeric(args[12])
human_reference_sample_file=as.character(args[13])
human_reference_legend_file=as.character(args[14])
human_reference_haplotype_file=as.character(args[15])
sampleName=as.character(args[16])

tempdir=outputdir
setwd(outputdir)
library("STITCH",lib.loc="./software/Rpackages-3.5.1/")
STITCH(
  bamlist = bamlist,
  reference = ref,
  outputdir = outputdir,
  method = "diploid",
  regenerateInput = TRUE,
  regionStart = regionStart,
  regionEnd = regionEnd,
  buffer = buffer,
  niterations = niterations,
  chr = chr,
  sampleNames_file = sampleName,
  inputBundleBlockSize = 100,
  reference_populations = c("CHB", "CHS", "CDX"),
  reference_haplotype_file = human_reference_haplotype_file,
  reference_sample_file = human_reference_sample_file,
  reference_legend_file = human_reference_legend_file,
  posfile = human_posfile, K = human_K, tempdir = tempdir, nCores = nCores, nGen = human_nGen)

```
usage: run_STITCH.sh
```bash

mkdir ./chr1_1_5000000

time ./software/01.Usr/bin/Rscript \
toolset/GWAS_creatCMD/STITCH.R \
./chr1_1_5000000 \
./bam_forge.list \
./toolset/reference/hg38/Homo_sapiens_assembly38.fasta \
./toolset/reference/1kg.easaf0.01/chr1.pos.txt \
10 16000 4 40 \
1 5000000 chr1 \
250000 \
./toolset/reference/liftover_easaf0.01/1000GP_Phase3.sample \
./toolset/reference/liftover_easaf0.01/chr1.legend.gz \
./toolset/reference/liftover_easaf0.01/chr1.hap.gz \
sampleid_forge.list

rm -rv ./chr1_1_5000000/input \
&& touch ./chr1_1_5000000/input.removed
```

### 1.3. Variants detection
usage: run_basevar.sh
```bash
BaseVarC basetype \
--input bam_forge.list \
--reference Homo_sapiens_assembly38.fasta \
--region chr12:1000001-1500000 \
--output ./chr12-1000001-1500000.10 \
--batch 10 --rerun 
```

### 1.4. Load in Plink & PCA
usage: run_plink_load.sh
```bash
plink2
--const-fid 0 \
--freq \
--hardy \
--make-pgen vzs \
--out WholeGenome2 \
--pca 5 \
--set-all-var-ids chr@:# \
--vcf WG.vcf dosage=DS
```

## 2. Statistical analyse
### 2.1. Genome wide association study (GWAS)
usage: run_GWAS.sh
```bash
./software/plink2 \
--pfile ./plink-format/STITCH/WholeGenome2 'vzs' \
--read-freq ./plink-format/STITCH/WholeGenome2.afreq \
--glm \
--maf 0.05 \
--hwe 1e-6 \
--geno 0.1 dosage \
--pheno ./phenotype_genaral/pheno.csv \
--pheno-quantile-normalize \
--covar ./plink-format/STITCH/basevar_PCA_age.eigenvec \
--covar-variance-standardize \
--out ./output/GENERAL/GWAS
```
### 2.2. Heritability and genetic correlation
code: h2.sh
```bash
glm_add=$1
name=$2
ofile=$3

# mkdir
mkdir ${ofile}/${name}
ofiles=${ofile}/${name}

# glm 2 ss
/share/app/python/3.8.6/bin/python ./toolset/glm_2_ss.py $1 ${ofiles}/${name}.ss

# ss 2 ldsc
/share/app/python/3.8.6/bin/python ./toolset/ss_2_ldsc.py ${ofiles}/${name}.ss ${ofiles}/${name}.ldsc

# get_size
sample_size=$(awk 'NR==2{print $2}' ${ofiles}/${name}.ldsc)
echo "${name}, sample_size: ${sample_size}"

# ldsc 2 munge
./envs/python27/bin/python2.7 \
./software/ldsc-master/munge_sumstats.py \
--sumstats ${ofiles}/${name}.ldsc \
--N $sample_size \
--out ${ofiles}/${name} \
--merge-alleles ./toolset/reference/wuhan_ss.snplist

# munge 2 h2
./envs/python27/bin/python2.7 \
./software/ldsc-master/ldsc.py \
--h2 ${ofiles}/${name}.sumstats.gz \
--ref-ld-chr ./toolset/LDSC/eas_ldscores/ \
--w-ld-chr ./toolset/LDSC/eas_ldscores/ \
--out ${ofiles}/${name}.h2
```
code: rg.sh
```bash
ss1=$1
ss2=$2
ofile=$3

./envs/python27/bin/python2.7 \
./software/ldsc-master/ldsc.py \
--rg $ss1,$ss2 \
--ref-ld-chr ./LDSC/eas_ldscores/ \
--w-ld-chr ./LDSC/eas_ldscores/ \
--out $ofile
```
usage:
```bash
bash glm_2_h2.sh ./output/GENERAL/GWAS.phenotype.add phenotype ./output/GENERAL/h2
bash rg.sh ./output/GENERAL/munge/phenotype1.ldsc.sumstats.gz ./output/GENERAL/munge/phenotype2.ldsc.sumstats.gz ./output/GENERAL/rg/RG_phenotype1_phenotype2
```

### 2.3. Gene expression analysis 
code: gtex.R
```R
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
```


### 2.4. Pathway enrichment analysis
code: PASCAL.sh
```bash
cd ./PASCAL/
./PASCAL/Pascal \
--pval=${ifile} \
--runpathway=on \
--custom=ASN \
--customdir=./PASCAL/resources/ASN/CHR \
--outsuffix=${ofile}
```


### 2.5. Partitioning heritability analysis
code: ph2.sh
```bash
./envs/python27/bin/python2.7 \
./software/ldsc-master/ldsc.py \
--h2-cts ${ofiles}/${name}.sumstats.gz \
--ref-ld-chr ./ldsc_reference/1000G_Phase3_EAS_baselineLD_v2.2_ldscore/baselineLD. \
--w-ld-chr ./ldsc_reference/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC. \
--frqfile-chr ./ldsc_reference/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC. \
--ref-ld-chr-cts ./ldsc_reference/Multi_tissue_gene_expr.EAS.ldcts \
--out ${ofiles}/${name}.ph2
```

### 2.6. Mendelian randomization analysis
code: MR.R
```R
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
```