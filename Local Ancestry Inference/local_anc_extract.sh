#!/bin/bash
OUT=`pwd`
SAMPLE=$1

RDATA=`ls ${OUT}/MOSAIC_RESULTS/${SAMPLE}*.RData | xargs -n 1 basename`

Rscript ./MOSAIC_anc_extract.R "${OUT}/MOSAIC_RESULTS/" $SAMPLE "${OUT}/MOSAIC_RESULTS/localanc/" $RDATA

for chr in {1..38}; do awk -v val="$chr" 'NR==1 {print; next} {$2=val; print}' "${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_${chr}_anc" | sed 's/ /\t/g' > "${OUT}/MOSAIC_RESULTS/localanc/temp_${SAMPLE}_${chr}_anc" && mv "${OUT}/MOSAIC_RESULTS/localanc/temp_${SAMPLE}_${chr}_anc" "${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_${chr}_anc"; done

cat ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_*_anc | awk '$5 == "European" && $7 == "European"' | awk '{print "chr"$2"_"$3}' > ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_European
cat ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_*_anc | awk '($5 == "European" && $7 == "PreContactDingo") || ($5 == "PreContactDingo" && $7 == "European")' | awk '{print "chr"$2"_"$3}' > ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_European_PreContactDingo
cat ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_*_anc | awk '$5 == "PreContactDingo" && $7 == "PreContactDingo"' | awk '{print "chr"$2"_"$3}' > ${OUT}/MOSAIC_RESULTS/localanc/${SAMPLE}_PreContactDingo