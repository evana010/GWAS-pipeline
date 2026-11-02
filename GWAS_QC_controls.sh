######## GWAS QC ###########
# load modules
module load PLINK/1.9b7.2-x86_64
module load R

# assign arguments
location=$1
filename=$2
startlocation=$3
strandfile=$4


################## Step 1. Strand Correction Wrayner ######################

# correct the strand with illumina strand correction avaialble on https://www.well.ox.ac.uk/~wrayner/strand/sourceStrand/index.html

#version_2/GSAMD-24v2-0_20024620_A1-b37.Source.strand
#version_3/GSAMD-24v3-0-EA_20034606_A1-b37.strand

bash /slade/home/sd718/All_cohort/eT1D/Population_startification/update_build.sh ${startlocation} ${strandfile} ${location}${filename}_wrayner

# add sex to .fam file
awk 'NR==FNR { lookup[$1] = $2; next } { if ($2 in lookup) $5 = lookup[$2] } { printf "%s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6 }'  /slade/projects/EXTEND/pituitary/analysis/QC/Controls/extend_pheno.txt ${location}${filename}_wrayner.fam > ${location}${filename}_sexupdated.fam
mv ${location}${filename}_sexupdated.fam ${location}${filename}_wrayner.fam

############################# STEP 2. Missingness #########################

# generate missingness report
plink --missing --bfile ${location}${filename}_wrayner --out ${location}${filename}
# output: plink.imiss and plink.lmiss, these files show respectively the missing SNPs per individual and the missing individuals per SNP

# plot missingness
Rscript --no-save plot_missing.R $location $filename
mv histimiss.pdf ${location}histimiss.pdf
mv histlmiss.pdf ${location}histlmiss.pdf


############ SNPs ############

# remove SNPs with 20% missingness
plink --bfile ${location}${filename}_wrayner --geno 0.2 --make-bed --out ${location}${filename}_geno20

# remove SNPs with 10% missing
plink --bfile ${location}${filename}_geno20 --geno 0.1 --make-bed --out ${location}${filename}_geno10

# remove SNPs with 5% missing
plink --bfile ${location}${filename}_geno10 --geno 0.05 --make-bed --out ${location}${filename}_geno5

# remove SNPs with 2% missingness (98% pass rate)
plink --bfile ${location}${filename}_geno5 --geno 0.02 --make-bed --out ${location}${filename}_geno2

########### Samples ###########

# remove samples with 20% missingness
plink --bfile ${location}${filename}_geno2 --mind 0.2 --make-bed --out ${location}${filename}_mind20

# remove samples with 10% missing
plink --bfile ${location}${filename}_mind20 --mind 0.1 --make-bed --out ${location}${filename}_mind10

# remove samles with 5% missing
plink --bfile ${location}${filename}_mind10 --mind 0.05 --make-bed --out ${location}${filename}_mind5

# remove samples with 2% missing
plink --bfile ${location}${filename}_mind5 --mind 0.02 --make-bed --out ${location}${filename}_mind2


############################ STEP 3. MAF ###################################

# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' ${startlocation}.bim > ${location}snp_1_22.txt
plink --bfile ${location}${filename}_mind2 --extract ${location}snp_1_22.txt --make-bed --out ${location}${filename}_1_22

# generate frequency report 
plink --bfile ${location}${filename}_1_22 --freq --out ${location}MAF_check

# plot freq
Rscript --no-save MAF_check.R $location
mv MAF_distribution.pdf ${location}MAF_distribution.pdf

# remove SNPs with a low MAF in sample set
plink --bfile ${location}${filename}_1_22 --maf 0.01 --make-bed --out ${location}${filename}_maf

############################# STEP 4. HWE ##################################

# generate HWE report
plink --bfile ${location}${filename}_maf --hardy --out ${location}hardy

# plot HWE
Rscript --no-save hwe.R $location
mv histhwe.pdf ${location}histhwe.pdf

# remove SNPs which deviate from the HWE
plink --bfile ${location}${filename}_maf --hwe 1e-6 --hwe-all --make-bed --out ${location}${filename}_hwe

##################### Step 5. Check for Duplicate IDs ########################

# check for duplicate IDs and add _repeat to any duplicates
Rscript --no-save duplicate_check.R $location $filename

# remove duplicates tagged with _repeat
awk '$2 ~ /repeat/ { print $1, $2 }' ${location}${filename}_hwe.fam > ${location}duplicates.txt
plink --bfile ${location}${filename}_hwe --remove ${location}duplicates.txt --make-bed --out ${location}${filename}_duprem

##################### Step 6. Sex Discrepencies #############################

# prune file with all chromosomes
plink --bfile ${location}${filename}_mind2 --exclude inversion_hg37.txt --range --indep-pairwise 1500 150 0.2 --out ${location}indepSNP_allchrom

# sex check on all chromosomes
plink --bfile ${location}${filename}_mind2 --extract ${location}indepSNP_allchrom.prune.in --check-sex --out ${location}sex_check

# remove sex discrepencies
grep "PROBLEM" ${location}sex_check.sexcheck| awk '{print$1,$2}'> ${location}sex_discrepancy.txt
plink --bfile ${location}${filename}_duprem --remove ${location}sex_discrepancy.txt --make-bed --out ${location}${filename}_sex


######################### STEP 7. Heterozygosity rate #########################
# high levels of heterozygosity within an individual might be an indication of low sample quality whereas low levels of heterozygosity may be due to inbreeding.
# high inversion regions on chromosomes 17 and 8 are exluded from this analysis as they will have a high heterozygosity rate

# removes SNPs within the inversion region
plink --bfile ${location}${filename}_sex --exclude inversion_hg37.txt --range --indep-pairwise 1500 150 0.2 --out ${location}indepSNP

# create heterozygosity report (excluding the inversion region)
plink --bfile ${location}${filename}_sex --extract ${location}indepSNP.prune.in --het --out ${location}het_check

# plot heterozygosity
Rscript --no-save heterozygosity_check.R $location
mv heterozygosity.pdf ${location}heterozygosity.pdf

# find samples 3SD away from the mean 
Rscript --no-save heterozygosity_outliers_list.R $location
sed 's/"// g' ${location}heterozygosity_outliers_list.txt | awk '{print$1, $2}'> ${location}heterozygosity_outliers.txt

# remove heterozygosity rate outliers
plink --bfile ${location}${filename}_sex --remove ${location}heterozygosity_outliers.txt --make-bed --out ${location}${filename}_het


########################## Step 8. Remove Related Samples ######################

# Check for relationships between individuals with a pihat > 0.2.
# 0.00: No sharing of alleles by descent (i.e., unrelated individuals).
# 0.25: Typically represents a first-degree relative, such as half-siblings, grandparent-grandchild, or uncle-niece.
# 0.50: Represents full siblings or parent-offspring.
# 1.00: Identical twins or duplicates.

plink --bfile ${location}${filename}_het --extract ${location}indepSNP.prune.in --genome --min 0.1875 --out ${location}pihat_min_related

# identify related samples (remove later for population stratification)
Rscript related_remove.R $location $filename

# plot relationships
Rscript --no-save Relatedness.R $location
mv relatedness.pdf ${location}relatedness.pdf
mv relatedness_hist.pdf ${location}relatedness_hist.pdf

plink --bfile ${location}${filename}_het --extract ${location}indepSNP.prune.in --genome --min 0.95 --out ${location}pihat_min_iden

# identify identical samples 
Rscript identical_remove.R $location $filename

# remove identical samples 
plink --bfile ${location}${filename}_het --remove ${location}related_remove.txt --make-bed --out ${location}${filename}_idenrem
