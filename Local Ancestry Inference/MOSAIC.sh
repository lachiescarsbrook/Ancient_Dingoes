#############################
# GLOBAL VARIABLE/DIRECTORY #
#############################
OUT=`pwd`
PLINK=$1

mkdir data
mkdir data/run
mkdir data/pops
mkdir data/geno

#Download and format recombination map files (Auton et al. 2013: https://doi.org/10.1371/journal.pgen.1003984)
git clone https://github.com/auton1/dog_recomb
RATES="${OUT}/dog_recomb/canFam3.1/maps" 
gunzip $RATES/*
for i in `seq 1 38`; do awk 'NR > 1 { print $2, $3, $4 }' $RATES/mark4_cleaned_chr${i}.cf3.1.sorted.txt > $RATES/rates.${i}; done

#Subset and format VCF file
plink --bfile ${PLINK} --dog --keep ${OUT}/mosaic_keep.txt --cm-map ${RATES}/rates.@ --mac 1 --bp-space 500 --make-bed --out ${OUT}/data/run/${PLINK}_sub
plink --bfile ${OUT}/data/run/${PLINK}_sub --dog --update-ids mosaic_keep.txt --make-bed --out ${OUT}/data/run/${PLINK}_mosaic
plink --bfile ${OUT}/data/run/${PLINK}_mosaic --dog --double-id --recode vcf --out ${OUT}/data/run/${PLINK}_mosaic
bcftools view -Oz -o ${OUT}/data/run/${PLINK}_mosaic.vcf.gz ${OUT}/data/run/${PLINK}_mosaic.vcf
bcftools index ${OUT}/data/run/${PLINK}_mosaic.vcf.gz

####################
# CREATE SNP FILES #
####################

#Separate into individual chromosomes and convert to MOSAIC snpfiles
bcftools query -f '%CHROM\n' ${OUT}/data/run/${PLINK}_mosaic.vcf.gz | sort | uniq > ${OUT}/data/run/chromosomes.txt

while read chr; do
    echo "Processing chromosome: $chr"
    bcftools view -r "$chr" ${OUT}/data/run/${PLINK}_mosaic.vcf.gz -Oz -o "${OUT}/data/run/chr${chr}.vcf.gz"
done < ${OUT}/data/run/chromosomes.txt

chromosomes=`cat ${OUT}/data/run/chromosomes.txt`
for CHR in $chromosomes; do
echo "Processing chromosome: $CHR"
#Convert into PLINK format
plink --vcf ${OUT}/data/run/chr${CHR}.vcf.gz --double-id --cm-map ${RATES}/rates.@ --dog --make-bed --out ${OUT}/data/run/chr${CHR}
#Convert bim files into MOSAIC snpfile format
cat ${OUT}/data/run/chr${CHR}.bim | awk '{print "chr"$1"_"$4 "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > ${OUT}/data/run/snpfile.${CHR};
done

#####################
# CREATE GENO FILES #
#####################

cat ${OUT}/data/run/${PLINK}_mosaic.fam | cut -f 1 -d " " | sort | uniq > test_pops
popset=`cat test_pops`
for i in $popset; do cat ${OUT}/data/run/${PLINK}_mosaic.fam | grep ${i} | cut -f 1,2 -d " " | sed 's/ /_/g'> ${OUT}/data/pops/${i}; done
export PLINK="Phased_Imputed_allchr1-38_CMmap"
cat test_pops | xargs -L 1 -P 0 ./geno_call.sh

#!/bin/bash
#OUT=`pwd`
#POP=$1
#PLINK=$2
#echo "Started ${POP}"
#bcftools view -S ${OUT}/data/pops/${POP} -Oz ${OUT}/data/run/${PLINK}_mosaic.vcf.gz -o ${OUT}/data/pops/${POP}.vcf.gz
#bcftools index -f ${OUT}/data/pops/${POP}.vcf.gz
#for i in `seq 1 38`; do bcftools view -r "${i}" -Oz ${OUT}/data/pops/${POP}.vcf.gz -o ${OUT}/data/pops/${POP}_${i}.vcf.gz; done
#echo "DONE ${POP} subset"
#for i in `seq 1 38`; do bcftools convert --hapsample --vcf-ids ${OUT}/data/pops/${POP}_${i}.vcf.gz -o ${OUT}/data/geno/${POP}_${i}; zcat ${OUT}/data/geno/${POP}_${i}.hap.gz | cut -f 6- -d " "| sed 's/ //g' | sed 's/*//g' > ${OUT}/data/geno/${POP}genofile.${i}; done
#echo "DONE ${POP} convert"

#Cleanup 
rm ${OUT}/data/geno/*.hap.gz ${OUT}/data/geno/*.samples

######################
# CREATE RATES FILES #
######################
#Create rates files
for i in `seq 1 38`;
do
	sites=$(cat ${OUT}/data/run/chr${i}.bim | wc -l)
	pos=$(cut -f 4 ${OUT}/data/run/chr${i}.bim)
	recom=$(cut -f 3 ${OUT}/data/run/chr${i}.bim)

	echo ":sites:$sites" > ${OUT}/data/run/rates.${i}
	echo $pos >> ${OUT}/data/run/rates.${i}
	echo $recom >> ${OUT}/data/run/rates.${i}
	echo "finished ${i}"
done

##############
# RUN MOSAIC #
##############
cat ${OUT}/data/pops/European ${OUT}/data/pops/PreContactAustralia | sed 's/_/ /' > sample.names
cut -f 3,4 mosaic_keep.txt | grep -v -e "European" -e "PreContact" | cut -f 2 > samples_to_run.txt

cat samples_to_run.txt | xargs -L 1 -P 0 ./mosaic_run.sh

#!/bin/bash
#OUT=`pwd`
#MOSAIC=$1

#mkdir data/run_${MOSAIC}
#cp data/run/rates.* data/run_${MOSAIC}/
#cp data/run/snpfile.* data/run_${MOSAIC}/
#cp data/geno/PreContactAustraliagenofile.* data/run_${MOSAIC}/
#cp data/geno/Europeangenofile.* data/run_${MOSAIC}/

#cp ${OUT}/sample.names ${OUT}/data/run_${MOSAIC}/sample.names
#cp ${OUT}/data/geno/${MOSAIC}genofile.{1..38} ${OUT}/data/run_${MOSAIC}/
#echo -e "\n${MOSAIC}\t${MOSAIC}" >> ${OUT}/data/run_${MOSAIC}/sample.names

#Rscript ${OUT}/MOSAIC/mosaic.R \
#	--nophase \
#	--ancestries 2 \
#	--chromosomes 1:38  \
#	-k 50 \
#	--maxcores 1 \
#	"${MOSAIC}" "${OUT}/data/run_${MOSAIC}"


###################################
# EXTRACT SNP-WISE LOCAL ANCESTRY #
###################################

mkdir ${OUT}/MOSAIC_RESULTS/localanc
cp ${OUT}/data/run/snpfile* ${OUT}/MOSAIC_RESULTS/

cat samples_to_run.txt | xargs -L 1 -P 0 ./local_anc_extract.sh

#!/bin/bash
#OUT=`pwd`
#SAMPLE=$1

#RDATA=`ls ${OUT}/MOSAIC_RESULTS/${SAMPLE}*.RData | xargs -n 1 basename`

#Assigns ancestry to each SNP based on a 90% probability threshold
#Rscript ./MOSAIC_anc_extract.R "${OUT}/MOSAIC_RESULTS/" $SAMPLE "${OUT}/MOSAIC_RESULTS/localanc/" $RDATA

#for chr in {1..38}; do awk -v val="$chr" 'NR==1 {print; next} {$2=val; print}' "${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_${chr}_anc" | sed 's/ /\t/g' > "${OUT}/MOSAIC_RESULTS/localanc/temp_${SAMPLE}_${chr}_anc" && mv "${OUT}/MOSAIC_RESULTS/localanc/temp_${SAMPLE}_${chr}_anc" "${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_${chr}_anc"; done

#Creates lists of SNPs for dingo-dingo, dingo-European and European-European ancestry blocks
#cat ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_*_anc | awk '$5 == "European" && $7 == "European"' | awk '{print "chr"$2"_"$3}' > ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_European
#cat ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_*_anc | awk '($5 == "European" && $7 == "PreContactDingo") || ($5 == "PreContactDingo" && $7 == "European")' | awk '{print "chr"$2"_"$3}' > ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_European_PreContactDingo
#cat ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_*_anc | awk '$5 == "PreContactDingo" && $7 == "PreContactDingo"' | awk '{print "chr"$2"_"$3}' > ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_PreContactDingo
