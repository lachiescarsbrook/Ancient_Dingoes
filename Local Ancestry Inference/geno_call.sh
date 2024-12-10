#!/bin/bash
OUT=`pwd`
POP=$1
PLINK=$2
echo "Started ${POP}"
bcftools view -S ${OUT}/data/pops/${POP} -Oz ${OUT}/data/run/${PLINK}_mosaic.vcf.gz -o ${OUT}/data/pops/${POP}.vcf.gz
bcftools index -f ${OUT}/data/pops/${POP}.vcf.gz
for i in `seq 1 38`; do bcftools view -r "${i}" -Oz ${OUT}/data/pops/${POP}.vcf.gz -o ${OUT}/data/pops/${POP}_${i}.vcf.gz; done
echo "DONE ${POP} subset"
for i in `seq 1 38`; do bcftools convert --hapsample --vcf-ids ${OUT}/data/pops/${POP}_${i}.vcf.gz -o ${OUT}/data/geno/${POP}_${i}; zcat ${OUT}/data/geno/${POP}_${i}.hap.gz | cut -f 6- -d " "| sed 's/ //g' | sed 's/*//g' > ${OUT}/data/geno/${POP}genofile.${i}; done
echo "DONE ${POP} convert"