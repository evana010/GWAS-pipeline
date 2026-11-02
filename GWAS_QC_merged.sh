######## GWAS QC ###########
# load modules
module load PLINK/1.9b7.2-x86_64
module load R

# assign arguments
location=$1
filename=$2
startlocation=$3


############################# STEP 2. Missingness #########################

# generate missingness report
plink --missing --bfile ${startlocation} --out ${location}${filename}
# output: plink.imiss and plink.lmiss, these files show respectively the missing SNPs per individual and the missing individuals per SNP

# plot missingness
Rscript --no-save plot_missing.R $location $filename
mv histimiss.pdf ${location}histimiss.pdf
mv histlmiss.pdf ${location}histlmiss.pdf


############ SNPs ############

# remove SNPs with 20% missingness
plink --bfile ${startlocation} --geno 0.2 --keep-allele-order --make-bed --out ${location}${filename}_geno20

# remove SNPs with 10% missing
plink --bfile ${location}${filename}_geno20 --geno 0.1 --keep-allele-order --make-bed --out ${location}${filename}_geno10

# remove SNPs with 5% missing
plink --bfile ${location}${filename}_geno10 --geno 0.05 --keep-allele-order --make-bed --out ${location}${filename}_geno5

# remove SNPs with 2% missingness (98% pass rate)
plink --bfile ${location}${filename}_geno5 --geno 0.02 --keep-allele-order --make-bed --out ${location}${filename}_geno2

########### Samples ###########

# remove samples with 20% missingness
plink --bfile ${location}${filename}_geno2 --mind 0.2 --keep-allele-order --make-bed --out ${location}${filename}_mind20

# remove samples with 10% missing
plink --bfile ${location}${filename}_mind20 --mind 0.1 --keep-allele-order --make-bed --out ${location}${filename}_mind10

# remove samles with 5% missing
plink --bfile ${location}${filename}_mind10 --mind 0.05 --keep-allele-order --make-bed --out ${location}${filename}_mind5

# remove samples with 2% missing
plink --bfile ${location}${filename}_mind5 --mind 0.02 --keep-allele-order --make-bed --out ${location}${filename}_mind2


############################ STEP 3. MAF ###################################

# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' ${startlocation}.bim > ${location}snp_1_22.txt
plink --bfile ${location}${filename}_mind2 --extract ${location}snp_1_22.txt --keep-allele-order --make-bed --out ${location}${filename}_1_22

# generate frequency report 
plink --bfile ${location}${filename}_1_22 --freq --out ${location}MAF_check

# plot freq
Rscript --no-save MAF_check.R $location
mv MAF_distribution.pdf ${location}MAF_distribution.pdf

# remove SNPs with a low MAF in sample set
plink --bfile ${location}${filename}_1_22 --maf 0.01 --keep-allele-order --make-bed --out ${location}${filename}_maf

############################# STEP 4. HWE ##################################

# generate HWE report
plink --bfile ${location}${filename}_maf --hardy --out ${location}hardy

# plot HWE
Rscript --no-save hwe.R $location
mv histhwe.pdf ${location}histhwe.pdf

# remove SNPs which deviate from the HWE
plink --bfile ${location}${filename}_maf --hwe 1e-6 --hwe-all --keep-allele-order --make-bed --out ${location}${filename}_hwe


######################### STEP 7. Heterozygosity rate #########################
# high levels of heterozygosity within an individual might be an indication of low sample quality whereas low levels of heterozygosity may be due to inbreeding.
# high inversion regions on chromosomes 17 and 8 are exluded from this analysis as they will have a high heterozygosity rate

# removes SNPs within the inversion region
plink --bfile ${location}${filename}_hwe --exclude inversion_hg37.txt --range --indep-pairwise 1500 150 0.2 --out ${location}indepSNP

# create heterozygosity report (excluding the inversion region)
plink --bfile ${location}${filename}_hwe --extract ${location}indepSNP.prune.in --het --out ${location}het_check

# plot heterozygosity
Rscript --no-save heterozygosity_check.R $location
mv heterozygosity.pdf ${location}heterozygosity.pdf

# find samples 3SD away from the mean 
Rscript --no-save heterozygosity_outliers_list.R $location
sed 's/"// g' ${location}heterozygosity_outliers_list.txt | awk '{print$1, $2}'> ${location}heterozygosity_outliers.txt

# remove heterozygosity rate outliers
plink --bfile ${location}${filename}_hwe --remove ${location}heterozygosity_outliers.txt --keep-allele-order --make-bed --out ${location}${filename}_het


########################## Step 8. Remove Related Samples ######################

# Check for relationships between individuals with a pihat > 0.2.
# 0.00: No sharing of alleles by descent (i.e., unrelated individuals).
# 0.25: Typically represents a first-degree relative, such as half-siblings, grandparent-grandchild, or uncle-niece.
# 0.50: Represents full siblings or parent-offspring.
# 1.00: Identical twins or duplicates.

########## Remove Related Samples ###########

plink --bfile ${location}${filename}_het --extract ${location}indepSNP.prune.in --genome --min 0.1875 --out ${location}pihat_min

# plot relationships
Rscript --no-save Relatedness.R $location
mv relatedness.pdf ${location}relatedness.pdf
mv relatedness_hist.pdf ${location}relatedness_hist.pdf

# remove samples 
Rscript relatedness_remove.R $location $filename
plink --bfile ${location}${filename}_het --remove ${location}related_remove.txt --keep-allele-order --make-bed --out ${location}${filename}_relatedrem


########## Remove Duplicate Samples #########

plink --bfile ${location}${filename}_het --extract ${location}indepSNP.prune.in --genome --min 0.95 --out ${location}pihat_min

# remove samples 
Rscript relatedness_remove.R $location $filename
plink --bfile ${location}${filename}_het --remove ${location}related_remove.txt --keep-allele-order --make-bed --out ${location}${filename}_idenrem


################ Flip SNPs ###################
plink --bfile ${location}${filename}_idenrem --a2-allele ${location}/snps_to_flip_a2.txt --make-bed --out ${location}${filename}_flipped


##################### Step 9. Convert to VCF for imputation ####################
module load BCFtools

# Output VCF 
plink --bfile ${location}${filename}_flipped --keep-allele-order --recode vcf --out ${location}${filename}
bgzip ${location}${filename}.vcf
bcftools index ${location}${filename}.vcf.gz


# split by chromosome 
for chr in $(bcftools query -f '%CHROM\n' EXTEND_merged_QC/CONTROLS_merged.vcf.gz | sort -u); do
    # Output file for the specific chromosome
    bcftools view -r $chr EXTEND_merged_QCC/CONTROLS_merged.vcf.gz -Oz -o EXTEND_merged_QC/vcf_by_chrom/EXTEND_${chr}.vcf.gz
done