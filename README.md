# Ultra-Low-Depth Sequencing Genetic Analysis

## Overview
This repository serves as a dedicated space for housing the codebase used in our publication for genetic analysis of ultra-low depth sequencing data. It is intended to provide researchers with access to the methodologies and algorithms employed in our study, facilitating further research and analysis in the field of genetic.

## License
The code within this repository is licensed under the [MIT License](./LICENSE.md). Please refer to the license file for more information on the terms and conditions of using and contributing to this project.

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
bash glm_2_h2.sh ./output/GENERAL/GWAS.GDM_AUC.glm.linear.add GDM_AUC ./output/GENERAL/h2
bash rg.sh ./output/GENERAL/munge/GDM_AUC.ldsc.sumstats.gz ./output/GENERAL/munge/pheno.ldsc.sumstats.gz ./output/GENERAL/rg/RG_GDM_AUC_pheno
```

### 2.3.MendelianRandomization(MR)
code: run_MR.R
```R
library(TwoSampleMR)
library(dplyr)
library(ieugwasr)

mr_results_list <- list()
result_list <- list()
r<-data.frame()

for (ef in exposure_files) {
  # Read exposure data
  exp_data <- read_exposure_data(ef, sep = "\t")

  exp_data <- clump_data(exp_data, clump_kb = 10000, clump_r2 = 0.2, clump_p1 = 1, clump_p2 = 1)
  exp_data<-exp_data %>% filter(SNP %in% clump_data$rsid)

  for (of in outcome_files) {
    # Read outcome data
    out_data <- read_outcome_data(of, sep = "\t")
  
    h <- harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)
    h <- h %>%
      mutate(
        beta.exposure = if_else(effect_allele.exposure != effect_allele.outcome, -beta.exposure, beta.exposure),
        effect_allele.exposure = if_else(effect_allele.exposure != effect_allele.outcome, other_allele.exposure, effect_allele.exposure),
        other_allele.exposure = if_else(effect_allele.exposure != effect_allele.outcome, effect_allele.exposure, other_allele.exposure)
      )
    
    dat <- mr(h, parameters = default_parameters(), method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

    exposure <- dat$exposure
    outcome <- dat$outcome
    method <- dat$method
    nsnp <- dat$nsnp
    b <- dat$b
    se <- dat$se
    pval <- dat$pval
    result_df <- data.frame(exposure, outcome, method, nsnp, b, se, pval)
    r<-rbind(r,result_df)
    result_list[[basename_prefix]] <- result_df

    write.table(
      x = result_df,
      file = output_file,
      sep = ",", 
      col.names = FALSE, 
      row.names = FALSE, 
      append = TRUE 
    )
  }
}
```


### 2.4.Transcriptome-wide association study(TWAS)
code: run_TWAS.R
```R
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
```

### 2.5.Drug target analysis
code: run_drug_target.sh
```bash
#step1_magma_annotation
magma \
--annotate window=35,10 \
--snp-loc 1kgp_eas.bim \
--gene-loc NCBI37.3.gene.loc \
--out Annotation

#step2_gene_analysis
#use=SNP column, P column ncol=N column
magma \
--bfile 1kgp_eas \
--pval GWAS_pheno.glm use=4,14 ncol=10 \
--gene-annot Annotation.genes.annot \
--gene-model multi \
--out GWAS_pheno

#step3_magma_Gene-set_analysis
magma \
--gene-results GWAS_pheno.genes.raw \
--set-annot drug_gene_set col=1,2 \
--out GWAS_pheno.drug

```
