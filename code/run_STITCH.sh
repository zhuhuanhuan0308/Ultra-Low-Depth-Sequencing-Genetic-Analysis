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