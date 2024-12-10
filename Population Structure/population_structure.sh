#All screening statistics (including sex determination) were calculated using the snakemake workflow CanID
#CanID also generated mitochondrial bam files, which were used in constructing phylogenies (see below)
git clone https://github.com/lachiescarsbrook/CanID.git #Follow the instructions

#Install required programs
conda create -n popstruct admixtools=7.0.2 eigensoft=8.0.0 plink=1.90b6.21 r booster=0.1.2 samtools=1.16 mafft iqtree=2.1.4
conda activate popstruct

#################
# VCF FILTERING	#
#################
#These PLINK files can be downloaded from the repository. DOI: 10.6084/m9.figshare.27256647
FILE="Pseudohaploid_allchr1-38"

#Download and format recombination map files (Auton et al. 2013: https://doi.org/10.1371/journal.pgen.1003984)
git clone https://github.com/auton1/dog_recomb
RATES="${OUT}/dog_recomb/canFam3.1/maps" 
gunzip $RATES/*
for i in `seq 1 38`; do awk 'NR > 1 { print $2, $3, $4 }' $RATES/mark4_cleaned_chr${i}.cf3.1.sorted.txt > $RATES/rates.${i}; done

#Filter vcf to exclude rare variants and sites in strong LD
plink --bfile ${FILE} --maf 0.01 --indep-pairwise 500 kb 1 0.5 --cm-map $RATES/rates.@ --dog --make-bed --out ${FILE}_maf0.01_LDpruned
plink --bfile ${FILE}_maf0.01_LDpruned --dog --recode --out ${FILE}_maf0.01_LDpruned

#########################
# EIGENSTRAT CONVERSION #
#########################

PLINK="${FILE}_maf0.01_LDpruned"
#Format 6th column to ensure compatibility with smartpca
awk '$6="1"' ${PLINK}.ped > ${PLINK}_mod.ped

echo genotypename: ${PLINK}_mod.ped > ${PLINK}.par.PED.EIGENSTRAT
echo snpname: ${PLINK}.map >> ${PLINK}.par.PED.EIGENSTRAT
echo indivname: ${PLINK}_mod.ped >> ${PLINK}.par.PED.EIGENSTRAT
echo outputformat: EIGENSTRAT >> ${PLINK}.par.PED.EIGENSTRAT
echo genotypeoutname: ${PLINK}.geno >> ${PLINK}.par.PED.EIGENSTRAT
echo snpoutname: ${PLINK}.snp >> ${PLINK}.par.PED.EIGENSTRAT
echo indivoutname: ${PLINK}.ind >> ${PLINK}.par.PED.EIGENSTRAT
echo familynames: NO >> ${PLINK}.par.PED.EIGENSTRAT

convertf -p ${PLINK}.par.PED.EIGENSTRAT

#Create new individual file with sample names in third column
awk '{print $1,$2,$1}' ${PLINK}.ind  > ${PLINK}.new.ind

##############################################################################################################

############
# SMARTPCA #
############
RUN="all_samples"
#Create new ind file containing population IDs
awk 'BEGIN { FS=OFS=" " } NR==FNR { a[$2]=$1; next } $1 in a { $3=a[$1] } 1' ${PLINK}.fam ${PLINK}.ind > ${PLINK}.group.ind
#Create a list of modern populations to be used in eigenvector construction (excluding outgroups)
less ${PLINK}.fam | cut -f 1 -d " " | grep -v "Ancient" | grep -v "jackal" | grep -v "Coyote01" | uniq > ${PLINK}.poplist.noancient
#Setup parfile
echo genotypename:     	${PLINK}.geno > par.PCA_${RUN}
echo snpname:           ${PLINK}.snp >> par.PCA_${RUN}
echo indivname:         ${PLINK}.group.ind >> par.PCA_${RUN}
echo evecoutname:       ${PLINK}.evec >> par.PCA_${RUN}
echo evaloutname:       ${PLINK}.eval >> par.PCA_${RUN}
echo poplistname:       ${PLINK}.poplist.noancient >> par.PCA_${RUN}
echo outliermode: 2 >> par.PCA_${RUN}
echo lsqproject:         YES >> $par.PCA_${RUN}
#Run
smartpca -p par.PCA_${RUN} > smartpca_${RUN}.out 

##############################################################################################################

###############
# Outgroup-f3 #
###############
#Download struct-f4 (https://bitbucket.org/plibradosanz/structf4/src/master/)
git clone https://bitbucket.org/plibradosanz/structf4.git #Be sure to install all neccesary dependancies
mkdir output_blocks

FILE="Pseudohaploid_allchr1-38"

#Convert into custom tped format
plink --bfile ${FILE} --dog --recode transpose --out ${FILE}
awk '{print $2,$1,$3,$4,$5,$2}' ${FILE}.tfam > ${FILE}_custom.tfam

#Creates Treemix-like files
perl Tped2Structf4.pl ${FILE}.tped ${FILE}_custom.tfam 5000000 output_blocks/files 1
ls output_blocks/file* > calcf4.filelist

#Determine position of outgroup taxa (Coyote)
head -n 1 calcf4.filelist > test_calcf4.filelist
EXAMPLE=`cat test_calcf4.filelist`
ALL=`awk '{print NF-1; exit}' ${EXAMPLE}`
OUTGREP=`grep 'Coyote' ${EXAMPLE} | awk '{for(i=1;i<=NF;i++){if($i=="Coyote"){print i}}}'`
OUTVAL=$((OUTGREP-2))

#Run calc-f3 to calculate every possible outgroup-f3 comparison
./Calc-f3 -i calcf4.filelist -o $OUTVAL -n $ALL -t 40 > out_f3_allcomparisons

#Converts output into matrix
python outf3_matrix.py "out_f3_allcomparisons"

##############################################################################################################

###########
# NJ TREE #
###########
Rscript -e 'install.packages("ape", repos = "https://cran.r-project.org")'

NAME="Pseudohaploid_allchr1-38"
OUT=`pwd`
mkdir trees

plink --bfile ${OUT}/${NAME} --distance square 1-ibs --dog --make-bed --out ${OUT}/${NAME}
cut -f 2 ${OUT}/${NAME}.bim > ${OUT}/site_list.txt
cut -f 2 ${OUT}/${NAME}.mdist.id |  perl -pe 's/\n/\t/g' | awk 'BEGIN{FS=OFS="ID\t"}{print value OFS $0}' > ${OUT}/matrix_head #Converts list to a tab-delimted line, and adds an 'ID' column at start (to account for offset)
cut -f 2 ${OUT}/${NAME}.mdist.id > ${OUT}/matrix_labels
paste ${OUT}/matrix_labels ${OUT}/${NAME}.mdist | cat ${OUT}/matrix_head - > ${OUT}/${NAME}.d

Rscript ${OUT}/bionj.R ${OUT}/${NAME}.d ${OUT}/${NAME}.tree

#bionj.R
	#args <- commandArgs(trailingOnly=TRUE)
	#library(ape)
	#distm<-read.table(args[1], header = TRUE, row.names = 1)
	#mdist.dist <- as.dist(distm)
	#nj.tree <- nj(mdist.dist)
	#max.dist <- max(node.depth.edgelength(nj.tree))
	#root.node <- which(node.depth.edgelength(nj.tree) == max.dist)
	#rooted.tree <- root(nj.tree, root.node)
	#write.tree(rooted.tree,file=args[2])

#Generate 100 replicate trees
for i in {1..100};
do 
	#Select 50,000 or 100,000 random sites to use to generate each tree
	echo -e "\n Subsetting Sites for Tree ${i}";
	cat site_list.txt | shuf | head -1000 > trees/sites_rep_${i}.txt; 
	
	#Calculate IBS distances between each sample in reduced dataset, and removes sites file afterwards
	echo -e "\n Creating IBS Distance Matrix for Tree ${i}";
	plink --bfile ${OUT}/${NAME} --distance square 1-ibs --extract trees/sites_rep_${i}.txt --dog --allow-no-sex --make-bed --out trees/${NAME}_rep_${i};
	rm trees/sites_rep_${i}.txt;
	
	#Formatting PLINK output
	cut -f 2 trees/${NAME}_rep_${i}.mdist.id |  perl -pe 's/\n/\t/g' | awk 'BEGIN{FS=OFS="ID\t"}{print value OFS $0}' > trees/matrix_head; #Converts list to a tab-delimted line, and adds an 'ID' column at start (to account for offset)
	cut -f 2 trees/${NAME}_rep_${i}.mdist.id > trees/matrix_labels;
	paste trees/matrix_labels trees/${NAME}_rep_${i}.mdist | cat trees/matrix_head - > trees/${NAME}_rep_${i}.d;

	#Create a tree specifying input/output files
	echo -e "\n Creating Tree ${i}";
	Rscript bionj.R trees/${NAME}_rep_${i}.d trees/${NAME}_rep_${i}.tree;

	#Create consensus tree with bootstrap
	cat trees/${NAME}_rep_${i}.tree | head -1 >> ${OUT}/${NAME}_bootstrap.tree;
	echo -e "\n Finished Tree ${i}";
done

cat trees/${NAME}_rep*.tree > ${OUT}/all_replicate.tree
booster -i ${OUT}/${NAME}_bootstrap.tree -b ${OUT}/all_replicate.tree -o ${OUT}/${NAME}_BOOSTER_tree.nwk

##############################################################################################################

##############
# MTDNA TREE #
##############
#Mitochondrial bam files were generated as part of CanID above
#Calculate summary statistics for mtDNA
echo "Sample Depth_of_Coverage Breadth_of_Coverage" > ${OUT}/coverage_statistics.txt

INDIR="path/to/CanID/mtDNA/bams"
OUT=`pwd`
MTBAM=`ls ${INDIR}*_mtDNA.bam | xargs -n 1 basename | sed 's/_mtDNA.bam//g'`

for i in $MTBAM; do
#Print sample name
	sample=$(echo ${i})
#Calculate depth of coverage (mitochondrial genome)
	depth=$(samtools depth ${INDIR}/${i}_mtDNA.bam | awk '{sum+=$3}END{print sum}' | awk '{print ($1/16730)}')
#Calculate breadth of coverage (mitochondrial genome)
	breadth=$(samtools depth ${INDIR}/${i}_mtDNA.bam | wc -l | awk '{print ($1/16730)*100 "%"}')
	echo "$sample $depth $breadth" >> ${OUT}/coverage_statistics.txt
done &

#Call consensus fasta files, and replace headers with sample names
#Calls fasta
for i in $MTBAM; do /home/lachie/downloads/samtools-1.16.1/samtools consensus -f FASTA -a -d 3 -c 0.75 --min-MQ 30 ${INDIR}/${i}_mtDNA.bam -o ${OUT}/${i}_mtDNA.fasta; done &
for i in $MTBAM; do cat ${i}_mtDNA.fasta | sed "s/chrM/${i}/g" > ${i}_mtDNA_consensus.fasta; done

#Create maximum-likelihood phylogeny
#Align sequences
#The fasta file can be downloaded from the repository. DOI: 10.6084/m9.figshare.27256647
mafft --auto --thread 20 Ancient_Dingoes_226.fasta > aligned_Ancient_Dingoes_226.fasta
#Create a maximum-likelihood tree from 1000 bootstrap replicates (following model optimization), specifying the outgroup (-o); with 1000 replicates for SH approximate likelihood ratio test (-alrt)
iqtree -s aligned_Ancient_Dingoes_226.fasta -st DNA -o 'NC008093.1_Coyote_USA' -B 1000 -ntmax 20 -T AUTO -m TEST 
