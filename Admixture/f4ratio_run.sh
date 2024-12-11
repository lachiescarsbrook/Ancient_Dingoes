#!/bin/bash
PLINK="/gpfs/scratch/pn29bi/ra65vij/GermanShepherd/Historic_German_Shepherd_maf0.01_thin"	
RUN=$1
SUFFIX=`echo $RUN | sed 's/f4ratio_input_//g'`
#Create parfile for qp3pop run
echo genotypename: ${PLINK}.geno > par.F4ratio_${SUFFIX}
echo snpname: ${PLINK}_f4ratio.snp >> par.F4ratio_${SUFFIX}
echo indivname: ${PLINK}_dates_reference.ind >> par.F4ratio_${SUFFIX}
echo popfilename: f4ratio_input_${SUFFIX} >> par.F4ratio_${SUFFIX}
echo numchrom: 38 >> par.F4ratio_${SUFFIX}
echo allsnps: YES >> par.F4ratio_${SUFFIX}
echo blgsize: 0.1 >> par.F4ratio_${SUFFIX}
echo printsd: YES >> par.F4STAT_${SUFFIX}
qpF4ratio -p par.F4ratio_${SUFFIX} > par.F4ratio_${SUFFIX}.out