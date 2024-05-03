#!/usr/bin/env bash
#SBATCH --mem=10GB
#SBATCH --job-name="QC-PRS"
#SBATCH --output=2-QC_target.log
#SBATCH -t 4-00:00  # time limit: (D-HH:MM) 


# ===> CHANGE HERE =================================================================================

PROJECTDIR="/home/s.defina/EGG" # project folder
GENDIR="/home/s.defina/GENR-Parents/Imputed/1000G_PhaseIIIv5" # location of target genetic data

# ==================================================================================================

# Set working directory to the project folder
cd $PROJECTDIR

# Create subfolder to store results
OUTDIR="$PROJECTDIR/results"
mkdir -p $OUTDIR

# ==================================================================================================
# Call vcftools to produce the heterozigosity statistics per chromosome and read into R to combine
echo " "
echo "======================== CLEANING TARGET FILE(S) ========================"
echo " "
echo "1. HETEROZIGOSITY"
#
#for i in {1..22}; do
#
#vcftools \
#--gzvcf $GENDIR/GENR_parents_Chr${i}.dose.vcf.gz \
#--het \
#--out $PROJECTDIR/het/het_chr${i}
#
#done
#
echo "Combining het files..."
#
#Rscript 2.1-make_het_tab.R $PROJECTDIR/het/
#
## Perform QC steps using plink 
#
echo "Done."
echo " "
echo "2. TARGET QC USING plink 1.9"

#ETHN="eur"  # "noneur"

#for ethn in $ETHN; do

for i in {1..22}; do
echo "=============================== chr ${i} ==============================="

# Call plink to perform all QC steps, output the list of individuals that pass the initial QC.
# Note:
# - id_list_eur.txt  file was created in 0-pheno_prep.R
# - excl_het_tab.txt file was created in 2.1-make_het_tab.R

plink19 \
--make-bed \
--vcf $GENDIR/GENR_parents_Chr${i}.dose.vcf.gz \
--mind 0.05 \
--geno 0.02 \
--hwe 0.000001 \
--maf 0.01 \
--keep $PROJECTDIR/pheno/id_list_eur.txt \
--remove $PROJECTDIR/het/excl_het_tab.txt \
--out $OUTDIR/chr${i}

echo "===================================================================="

done

# Note: the (cleaned) dataset is split into chromosomes
echo " "
echo "Merging (cleaned) chromosomes..."
cd $OUTDIR
seq 1 22 | xargs -I _num echo chr_num > mergelist.txt
plink2 --pmerge-list mergelist.txt 'bfile' --make-bed --pmerge-output-vzs --out all_autosomal_eur

# Clean up 
rm chr*

#done