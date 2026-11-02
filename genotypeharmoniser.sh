
######################### Step 9. Genotype Harmoniser ################################

genotype_harmoniser() {
	# run GH to change SNP IDs to a reference 
	/slade/projects/EXTEND/pituitary/analysis/QC/GenotypeHarmonizer/GenotypeHarmonizer.sh --update-reference-allele --update-id \
		--input ${1}_idenrem --inputType PLINK_BED \
		--output ${1}_genotypeharmoniser_1kg --outputType PLINK_BED \
		--ref /slade/local/pituitary/analysis/Population_stratification/1kg_ref/1kg_phase1_all --refType PLINK_BED
}

genotype_harmoniser Batch0.5/GSA2016_250_025
genotype_harmoniser Batch1/GC-KM-10290-All-Sample-ST
genotype_harmoniser Controls/GSA2018_283_025/GSA2018_283_025
genotype_harmoniser Controls/GSA2018-377-025/GSA2018-377-025
genotype_harmoniser Controls/GSA2016_16_025/GSA2016_16_025
genotype_harmoniser Batch2/PLINK_060525_merged_coloradorem