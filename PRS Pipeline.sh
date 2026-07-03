#!/usr/bin/env bash

# Polygenic Risk Score (PRS) Pipeline

# Software:

#   - PLINK 1.9

#   - PLINK2

#   - PRSice-2

#   - PRS-CS

#   - PRS-CSx

#   - Python 3

#   - R


set -e

set -o pipefail

# Input files

PATH = /mnt/Storage/Kathir/

GENOTYPE=sczctrl

RSID_MAP=sczctrl_filtered_chr_bp_sorted_rsID.txt

PGC_SNPS=PGC3_SCZ_SNP_markers.txt

GWAS=PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv

PRSICE=PRSice_linux

PRSICE_SCRIPT=PRSice.R

# 1. Genotype QC

plink \

    --bfile ${GENOTYPE} \

    --geno 0.05 \

    --mind 0.02 \

    --maf 0.01 \

    --hwe 1e-6 \

    --make-bed \

    --out sczctrl_filtered


# 2. Add rsIDs

sh add_rsID_bim.sh \

    ${RSID_MAP} \

    sczctrl_filtered.bim \

    sczctrl_filtered_rsID.bim


# 3. Extract SNPs common with PGC

cut -f2 sczctrl_filtered_rsID.bim > our_snps.txt

sort our_snps.txt > our_snps.sorted

sort ${PGC_SNPS} > pgc_snps.sorted

comm -12 pgc_snps.sorted our_snps.sorted > common_snps.txt

plink \

    --bfile sczctrl_filtered \

    --extract common_snps.txt \

    --make-bed \

    --out pgc_common


# 4. Remove related individuals

plink \

    --bfile pgc_common \

    --genome \

    --out relatedness

awk '$10>0.2 {print $1,$2}' relatedness.genome > related_remove.txt

plink \

    --bfile pgc_common \

    --remove related_remove.txt \

    --make-bed \

    --out pgc_common_unrelated


# 5. LD pruning


plink \

    --bfile pgc_common_unrelated \

    --indep-pairwise 200 50 0.2 \

    --out pgc_common

plink \

    --bfile pgc_common_unrelated \

    --extract pgc_common.prune.in \

    --make-bed \

    --out pgc_common_LD


# 6. Principal Component Analysis

plink \

    --bfile pgc_common_LD \

    --pca 20 \

    --out PCA


# 7. Prepare GWAS summary statistics

awk -F'\t' '

BEGIN{OFS="\t"}

$1 !~ /^#/ {print}

' ${GWAS} > GWAS_noheader.tsv



awk -F'\t' '

BEGIN{OFS="\t"}

NR==1{

print "SNP","CHR","BP","A1","A2","BETA","P"

next

}

{

print $2,$1,$3,$4,$5,$9,$11

}

' GWAS_noheader.tsv > GWAS_sumstats.txt

# 8. Prepare PRS-CS summary statistics

awk '

BEGIN{OFS="\t"}

NR==1{

print "SNP","A1","A2","BETA","P"

next

}

{

print $1,$4,$5,$6,$7

}

' GWAS_sumstats.txt > GWAS_PRScs.txt

# 9. PRSice

Rscript ${PRSICE_SCRIPT} \

    --prsice ${PRSICE} \

    --base GWAS_sumstats.txt \

    --target pgc_common_unrelated \

    --pheno phenotype.txt \

    --cov PCA.txt \

    --base-maf MAF:0.01 \

    --stat BETA \

    --binary-target T \

    --pheno-col Phenotype \

    --cov-col PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \

    --clump-kb 250 \

    --clump-r2 0.1 \

    --clump-p 1 \

    --score std \

    --out prsice_output

# 10. PRS-CS

mkdir -p prscs_output

for chr in {1..22}

do

python PRScs.py \

    --ref_dir ldblk_1kg_eur \

    --bim_prefix pgc_common_unrelated \

    --sst_file GWAS_PRScs.txt \

    --n_gwas 170114 \

    --chrom ${chr} \

    --out_dir prscs_output

done

# 11. Merge PRS-CS weights

cat prscs_output/*chr*.txt > prscs_weights.txt

# 12. Calculate PRS

plink2 \

    --bfile pgc_common_unrelated \

    --score prscs_weights.txt 2 4 6 header cols=+scoresums ignore-dup-ids \

    --out prscs_scores


# 13. PRS-CSx

mkdir -p PRScsx

for chr in {1..22}

do

python PRScsx.py \

    --ref_dir ldblk_1kg \

    --bim_prefix pgc_common_unrelated \

    --sst_file \

PGC3_SCZ_wave3.european.autosome_clean.txt,PGC3_SCZ_wave3.asian.autosome.public.v3_clean.txt \

    --n_gwas 127691,30075 \

    --pop EUR,EAS \

    --chrom ${chr} \

    --phi 1e-2 \

    --out_dir PRScsx \

    --out_name ctrl_scz_PRScsx

done

echo "========================================"

echo "      PRS Pipeline Completed"

echo "========================================"
