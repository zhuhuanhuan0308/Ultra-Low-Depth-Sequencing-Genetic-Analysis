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