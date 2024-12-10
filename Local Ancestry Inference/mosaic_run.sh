#!/bin/bash
OUT=`pwd`
MOSAIC=$1

mkdir data/run_${MOSAIC}
cp data/run/rates.* data/run_${MOSAIC}/
cp data/run/snpfile.* data/run_${MOSAIC}/
cp data/geno/PreContactAustraliagenofile.* data/run_${MOSAIC}/
cp data/geno/Europeangenofile.* data/run_${MOSAIC}/

cp ${OUT}/sample.names ${OUT}/data/run_${MOSAIC}/sample.names
cp ${OUT}/data/geno/${MOSAIC}genofile.{1..38} ${OUT}/data/run_${MOSAIC}/
echo -e "\n${MOSAIC}\t${MOSAIC}" >> ${OUT}/data/run_${MOSAIC}/sample.names

Rscript ${OUT}/MOSAIC/mosaic.R \
	--nophase \
	--ancestries 2 \
	--chromosomes 1:38  \
	-k 50 \
	--maxcores 30 \
	"${MOSAIC}" "${OUT}/data/run_${MOSAIC}"