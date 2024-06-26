./envs/python27/bin/python2.7 \
./software/ldsc-master/ldsc.py \
--h2-cts ${ofiles}/${name}.sumstats.gz \
--ref-ld-chr ./ldsc_reference/1000G_Phase3_EAS_baselineLD_v2.2_ldscore/baselineLD. \
--w-ld-chr ./ldsc_reference/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC. \
--frqfile-chr ./ldsc_reference/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC. \
--ref-ld-chr-cts ./ldsc_reference/Multi_tissue_gene_expr.EAS.ldcts \
--out ${ofiles}/${name}.ph2