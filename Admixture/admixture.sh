

conda create -n admixture admixture=1.3.0 admixtools=7.0.2

#############
# ADMIXTURE #
#############
PLINK="${FILE}_maf0.01_LDpruned"

#Run ADMIXTURE on all dogs
cat ${PLINK}.fam | cut -f 1,2 -d " " | grep -v "Coyote" > exclude_coyote.txt
plink --bfile ${PLINK} --keep exclude_coyote.txt --mac 1 --dog --make-bed --out ADMIXTURE_nocoyote
less ADMIXTURE_nocoyote.fam | cut -f 2 -d " " > ADMIXTURE_nocoyote_samples.txt
for i in `seq 2 6`; do admixture --cv -B 50 ADMIXTURE_nocoyote.bed ${i} > ADMIXTURE_nocoyote_log_${i}.out; done

#Run ADMIXTURE on dingoes and European dogs only
cat ${FILE}.fam | cut -f 1,2 -d " " | grep -e "European" -e "Dingo" -e "NewGuinea" > dingoes_european.txt
plink --bfile ${PLINK} --keep dingoes_european.txt --mac 1 --dog --make-bed --out ADMIXTURE_dingoesEU
less ADMIXTURE_dingoesEU.fam | cut -f 2 -d " " > ADMIXTURE_dingoesEU_samples.txt
for i in `seq 2 4`; do admixture --cv -B 50 ADMIXTURE_dingoesEU.bed ${i} > ADMIXTURE_dingoesEU_log_${i}.out; done

##############################################################################################################

#########
# DSTAT #
#########
PLINK="${FILE}_maf0.01_LDpruned"

cat ${PLINK}.fam | grep -e "European" | grep -v "Historic" | cut -f 2 -d " " > european.txt
cat european.txt | xargs -L 1 -P 0 ./dstat_setup.sh

#Split into smaller chunks
split -l 200 qpDstat_input qpDstat_input_
ls qpDstat_input_* > Dstat_blocks

#Runs the outgroup F3 for each block
cat Dstat_blocks | xargs -L 1 -P 0 ./dstat_run.sh
wait

#!/bin/bash
#FILE="Pseudohaploid_allchr1-38"
#PLINK="${FILE}_maf0.01_LDpruned"
#RUN=$1
#SUFFIX=`echo $RUN | sed 's/qpDstat_input_//g'`
#Create qpDstat parfile
#echo genotypename: ${PLINK}.geno > par.F4STAT_${SUFFIX}
#echo snpname: ${PLINK}.snp >> par.F4STAT_${SUFFIX}
#echo indivname: ${PLINK}.new.ind >> par.F4STAT_${SUFFIX}
#echo popfilename: qpDstat_input_${SUFFIX} >> par.F4STAT_${SUFFIX}
#echo allsnps: NO >> par.F4STAT_${SUFFIX}
#echo numchrom: 38 >> par.F4STAT_${SUFFIX}
#echo printsd: YES >> par.F4STAT_${SUFFIX}
#qpDstat -p par.F4STAT_${SUFFIX} > par.F4STAT_${SUFFIX}.out

cat par.F4STAT_* | egrep 'result' | sed 's/  */\t/g' | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > all.results
       
############
# F4-ratio #
############

PLINK="${FILE}_maf0.01_LDpruned"
cat ${PLINK}.fam | grep -e "European" | grep -v "Historic" | cut -f 2 -d " " > european.txt
cat european.txt | xargs -L 1 -P 0 ./f4ratio_setup.sh

#Sets false physical positions and genetic distances to allow block jacknife for 100 SNPs
awk '{$4=10000*NR}1' ${PLINK}.snp | awk '{$3=0}1' > ${PLINK}_f4ratio.snp

split -l 200 f4ratio_input f4ratio_input_
ls f4ratio_input_* > f4ratio_blocks

#Runs the outgroup F3 for each block
cat f4ratio_blocks | xargs -L 1 -P 0 ./f4ratio_run.sh
wait

#!/bin/bash
#PLINK="/gpfs/scratch/pn29bi/ra65vij/GermanShepherd/Historic_German_Shepherd_maf0.01_thin"	
#RUN=$1
#SUFFIX=`echo $RUN | sed 's/f4ratio_input_//g'`
#Create parfile for qp3pop run
#echo genotypename: ${PLINK}.geno > par.F4ratio_${SUFFIX}
#echo snpname: ${PLINK}_f4ratio.snp >> par.F4ratio_${SUFFIX}
#echo indivname: ${PLINK}_dates_reference.ind >> par.F4ratio_${SUFFIX}
#echo popfilename: f4ratio_input_${SUFFIX} >> par.F4ratio_${SUFFIX}
#echo numchrom: 38 >> par.F4ratio_${SUFFIX}
#echo allsnps: YES >> par.F4ratio_${SUFFIX}
#echo blgsize: 0.1 >> par.F4ratio_${SUFFIX}
#echo printsd: YES >> par.F4STAT_${SUFFIX}
#qpF4ratio -p par.F4ratio_${SUFFIX} > par.F4ratio_${SUFFIX}.out

egrep 'result' par.F4ratio_*.out  | sed 's/  */\t/g' | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}'  > par.F4ratio.results