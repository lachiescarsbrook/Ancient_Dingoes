#!/bin/bash
FILE="Pseudohaploid_allchr1-38"
PLINK="${FILE}_maf0.01_LDpruned"
RUN=$1
SUFFIX=`echo $RUN | sed 's/qpDstat_input_//g'`
#Create qpDstat parfile
echo genotypename: ${PLINK}.geno > par.F4STAT_${SUFFIX}
echo snpname: ${PLINK}.snp >> par.F4STAT_${SUFFIX}
echo indivname: ${PLINK}.new.ind >> par.F4STAT_${SUFFIX}
echo popfilename: qpDstat_input_${SUFFIX} >> par.F4STAT_${SUFFIX}
echo allsnps: NO >> par.F4STAT_${SUFFIX}
echo numchrom: 38 >> par.F4STAT_${SUFFIX}
echo printsd: YES >> par.F4STAT_${SUFFIX}
qpDstat -p par.F4STAT_${SUFFIX} > par.F4STAT_${SUFFIX}.out