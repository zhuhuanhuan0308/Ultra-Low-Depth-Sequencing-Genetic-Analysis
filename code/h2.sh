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
