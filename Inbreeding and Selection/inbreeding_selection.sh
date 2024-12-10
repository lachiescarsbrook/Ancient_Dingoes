mamba create -n inbreedselection bedtools=2.30.0 plink=1.90b6.21 vcftools=0.1.16
mamba activate inbreedselection

####### 
# ROH #
#######
FILE="Phased_Imputed_allchr1-38_CMmap"

#Convert VCF into PLINK format
plink --bfile ${FILE} --hardy --dog --out ${FILE}

#Creates a file for summary statistics
awk 'END {print "Imputed_SNPs: " NR}' ${FILE}.bim > summary.txt
awk '{sum += $7} END {print "Imputed_Mean_Heterozygosity_Observed: " sum/NR}' ${FILE}.hwe >> summary.txt
awk '{sum += $8} END {print "Imputed_Mean_Heterozygosity_Expected: " sum/NR}' ${FILE}.hwe >> summary.txt

#Calculate L (Minimum SNPs per ROH) and t (threshold)
Ns=$(cat summary.txt | grep "Imputed_SNPs" | cut -f 2 -d " ")
ind=$(cat ${FILE}.fam | wc -l)
alpha=0.05
het=$(cat summary.txt | grep "Imputed_Mean_Heterozygosity_Observed" | cut -f 2 -d " ")
awk "{print (log($alpha/($Ns*$ind))/(log(1-$het)))}"

#Estimate ROH for each individual 
plink --bfile ${FILE} --homozyg --homozyg-window-snp 100 --homozyg-window-threshold 0.02 \
	--homozyg-gap 250 --homozyg-snp 50 --homozyg-kb 130 --homozyg-density 50 \
	--homozyg-window-het 3 --homozyg-het 2 --dog --make-bed --out ${FILE}_ROH

#All further calculations and plotting was performed using the RScript ROH Plot

#############################################################################################################

#############
# SELECTION #
#############
#Calculating FST
plink --bfile ${FILE} --dog --recode vcf-iid --out ${FILE}
#Windowed
vcftools --vcf ${FILE}.vcf --weir-fst-pop precontact.txt --weir-fst-pop postcontact.txt --fst-window-size 50000 --fst-window-step 10000 --out pre_post_dingoes
#SNP-wise
vcftools --vcf ${FILE}.vcf --weir-fst-pop precontact.txt --weir-fst-pop postcontact.txt --out pre_post_dingoes

#Determine genes under selection
#Download CanFam3.1 annotation
wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/genes/canFam3.ncbiRefSeq.gtf.gz
gunzip canFam3.ncbiRefSeq.gtf.gz
#Format file
cat canFam3.ncbiRefSeq.gtf | grep "gene" | awk '{print $1, $4,$5,$10}' | tr -d '";' | sed 's/ /\t/g' > canFam3.ncbiRefSeq.bed

#Run BEDTOOLS intersect
#region.bed is generated using the RScript ROH Plot
bedtools intersect -a region.bed -b canFam3.ncbiRefSeq.bed -wa -wb > snps_genes_overlap_ncbiRefSeq_wolfdogs.bed

