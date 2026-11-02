#!/usr/bin/env bash
conda activate myenv

# filter R squared
for i in {1..22}; do
    FILEPATH="/slade/projects/EXTEND/pituitary/analysis/Imputation/chr_${i}/chr${i}.dose.vcf.gz"
    OUTPUT_FILE="chr_${i}_filtered.vcf.gz"
    if [[ -f $FILEPATH ]]; then
        bcftools view -i 'INFO/R2>=0.3' $FILEPATH -O z -o $OUTPUT_FILE
    fi
done

# alternative method using multiple tmux
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_1/chr1.dose.vcf.gz -O z -o chr_1_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_2/chr2.dose.vcf.gz -O z -o chr_2_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_3/chr3.dose.vcf.gz -O z -o chr_3_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_4/chr4.dose.vcf.gz -O z -o chr_4_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_5/chr5.dose.vcf.gz -O z -o chr_5_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_6/chr6.dose.vcf.gz -O z -o chr_6_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_7/chr7.dose.vcf.gz -O z -o chr_7_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_8/chr8.dose.vcf.gz -O z -o chr_8_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_9/chr9.dose.vcf.gz -O z -o chr_9_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_10/chr10.dose.vcf.gz -O z -o chr_10_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_11/chr11.dose.vcf.gz -O z -o chr_11_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_12/chr12.dose.vcf.gz -O z -o chr_12_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_13/chr13.dose.vcf.gz -O z -o chr_13_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_14/chr14.dose.vcf.gz -O z -o chr_14_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_15/chr15.dose.vcf.gz -O z -o chr_15_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_16/chr16.dose.vcf.gz -O z -o chr_16_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_17/chr17.dose.vcf.gz -O z -o chr_17_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_18/chr18.dose.vcf.gz -O z -o chr_18_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_19/chr19.dose.vcf.gz -O z -o chr_19_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_20/chr20.dose.vcf.gz -O z -o chr_20_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_21/chr21.dose.vcf.gz -O z -o chr_21_filtered.vcf.gz
bcftools view -i 'INFO/R2>=0.3' /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_22/chr22.dose.vcf.gz -O z -o chr_22_filtered.vcf.gz


#convert filtered to plink
for i in {1..22}; do
    FILEPATH="/slade/projects/EXTEND/pituitary/analysis/Imputation_QC/chr_${i}_filtered.vcf.gz"
    OUTPUT_FILE="chr_${i}_filtered"
    if [[ -f $FILEPATH ]]; then
        plink --vcf $FILEPATH --make-bed --out $OUTPUT_FILE
    fi
done

# alternative method using multiple tmux
plink --vcf chr_1_filtered.vcf.gz --make-bed --out chr_1_filtered
plink --vcf chr_2_filtered.vcf.gz --make-bed --out chr_2_filtered
plink --vcf chr_3_filtered.vcf.gz --make-bed --out chr_3_filtered
plink --vcf chr_4_filtered.vcf.gz --make-bed --out chr_4_filtered
plink --vcf chr_5_filtered.vcf.gz --make-bed --out chr_5_filtered
plink --vcf chr_6_filtered.vcf.gz --make-bed --out chr_6_filtered
plink --vcf chr_7_filtered.vcf.gz --make-bed --out chr_7_filtered
plink --vcf chr_8_filtered.vcf.gz --make-bed --out chr_8_filtered
plink --vcf chr_9_filtered.vcf.gz --make-bed --out chr_9_filtered
plink --vcf chr_10_filtered.vcf.gz --make-bed --out chr_10_filtered
plink --vcf chr_11_filtered.vcf.gz --make-bed --out chr_11_filtered
plink --vcf chr_12_filtered.vcf.gz --make-bed --out chr_12_filtered
plink --vcf chr_13_filtered.vcf.gz --make-bed --out chr_13_filtered
plink --vcf chr_14_filtered.vcf.gz --make-bed --out chr_14_filtered
plink --vcf chr_15_filtered.vcf.gz --make-bed --out chr_15_filtered
plink --vcf chr_16_filtered.vcf.gz --make-bed --out chr_16_filtered
plink --vcf chr_17_filtered.vcf.gz --make-bed --out chr_17_filtered
plink --vcf chr_18_filtered.vcf.gz --make-bed --out chr_18_filtered
plink --vcf chr_19_filtered.vcf.gz --make-bed --out chr_19_filtered
plink --vcf chr_20_filtered.vcf.gz --make-bed --out chr_20_filtered
plink --vcf chr_21_filtered.vcf.gz --make-bed --out chr_21_filtered
plink --vcf chr_22_filtered.vcf.gz --make-bed --out chr_22_filtered


# extract R2 values for all snps
echo -e "CHROM\tPOS\tID\tREF\tALT\tR2" > all_SNPs_R2.txt
for i in {1..22}; do
    zcat /slade/projects/EXTEND/pituitary/analysis/Imputation/chr_${i}/chr${i}.info.gz | 
    awk 'BEGIN { OFS="\t"; FS="\t" }  
         NR > 1 {
             split($8, info_fields, ";");  
             for (j in info_fields) { 
                 split(info_fields[j], kv, "="); 
                 if (kv[1] == "R2") print $1, $2, $3, $4, $5, kv[2]; 
             } 
         }' >> all_SNPs_R2.txt
done


# rename snps with full stop
for file in *.bim; do
  awk '{if ($2 == ".") $2 = $1 ":" $4 ":" $5 ":" $6}1' OFS="\t" "$file" > temp_file && mv temp_file "$file"
done

# merge all (produce missnp)
plink --bfile chr_1_filtered --merge-list merge_list.txt --make-bed --out imputed_merged


# remove missnp
for file in *.bed; do
  # Extract base name without extension (e.g., data for data.bed, data.bim, data.fam)
  base=${file%.bed}
  
  # Remove SNPs listed in missnp from the binary files
  plink --bfile "$base" --exclude imputed_merged-merge.missnp --make-bed --out "${base}_missnprem"
done

# merge all
plink --bfile chr_1_filtered_missnprem --merge-list merge_list_2.txt --make-bed --out imputed_merged


# add sex and phenotype
cp /slade/projects/EXTEND/pituitary/analysis/QC/Cases_controls_merged/CASES_and_CONTROLS_merged_idenrem.fam imputed_merged.fam

# remove AIPneg unaffected cases
plink --bfile imputed_merged --remove unaffected_AIPneg_REMOVE.txt --make-bed --out imputed_merged_unaffectedrem

# QC imputed merged data
plink --bfile imputed_merged_unaffectedrem --geno 0.02 --make-bed --out imputed_merged_geno
plink --bfile imputed_merged_geno --mind 0.02 --make-bed --out imputed_merged_mind
plink --bfile imputed_merged_mind --maf 0.01 --make-bed --out imputed_merged_maf
plink --bfile imputed_merged_maf --hwe 1e-6 --make-bed --out imputed_merged_hwe


# check for relatedness
plink --bfile /slade/projects/EXTEND/pituitary/analysis/QC/Cases_controls_merged/CASES_and_CONTROLS_merged_idenrem --extract /slade/projects/EXTEND/pituitary/analysis/QC/Cases_controls_merged/indepSNP.prune.in --genome --min 0.1875 --out pihat_min_related
Rscript related_remove.R /slade/projects/EXTEND/pituitary/analysis/QC/Cases_controls_merged/CASES_and_CONTROLS_merged

