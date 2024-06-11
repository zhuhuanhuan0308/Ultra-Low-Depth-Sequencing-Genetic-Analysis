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