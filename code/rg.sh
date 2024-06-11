ss1=$1
ss2=$2
ofile=$3

./envs/python27/bin/python2.7 \
./software/ldsc-master/ldsc.py \
--rg $ss1,$ss2 \
--ref-ld-chr ./LDSC/eas_ldscores/ \
--w-ld-chr ./LDSC/eas_ldscores/ \
--out $ofile