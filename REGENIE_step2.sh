#!/bin/bash

module purge
module load GCC/9.3.0 

/slade/projects/EXTEND/pituitary/analysis/Association_testing/regenie \
  --step 2 \
  --bed "/slade/projects/EXTEND/pituitary/analysis/Association_testing/AIP_neg" \
  --covarFile "/slade/projects/EXTEND/pituitary/analysis/Association_testing/covariates_all.txt" \
  --phenoFile "/slade/projects/EXTEND/pituitary/analysis/Association_testing/pheno_all.txt"  \
  --bsize 200 \
  --bt \
  --firth --approx \
  --pThresh 0.01 \
  --pred /slade/projects/EXTEND/pituitary/analysis/Association_testing/GWAS_step1_output_allinc_pred.list \
  --out GWAS_step2_output_allinc