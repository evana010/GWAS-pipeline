#!/bin/bash

module purge
module load GCC/9.3.0 

/slade/projects/EXTEND/pituitary/analysis/Association_testing/regenie \
   --step 1 \
   --bed "/slade/projects/EXTEND/pituitary/analysis/Association_testing/AIP_neg_european_preimp" \
   --covarFile "/slade/projects/EXTEND/pituitary/analysis/Association_testing/covariates_all.txt" \
   --phenoFile "/slade/projects/EXTEND/pituitary/analysis/Association_testing/pheno_all.txt"  \
   --bsize 1000 \
   --threads 16 \
   --bt \
   --loocv \
   --maxCatLevels 10 \
   --catCovarList Sex \
   --covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
   --out GWAS_step1_output_allinc_EURONLY



