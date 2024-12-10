The scripts contained within this directory were used to perform local ancestry inference in MOSAIC, and then extract regions of ancestry for use in outgroup-f3 comparisons.

First, clone this repository, and install all of the required packages. Make sure your working directory is `Local Ancestry Inference`

# CONDA ENVIRONMENT #
conda create -n mosaic plink=1.90b6.21 bcftools=1.21 r
conda activate mosaic

# MOSAIC INSTALL #
git clone https://github.com/mststats/MOSAIC.git
wget https://maths.ucd.ie/~mst/MOSAIC/MOSAIC_1.5.1.tar.gz

Before starting, you must ensure all of the dependencies for MOSAIC have been installed, and the Phased_Imputed VCF is in the working directory. Then run:

bash ./MOSAIC.sh "Phased_Imputed_allchr1-38_CMmap"