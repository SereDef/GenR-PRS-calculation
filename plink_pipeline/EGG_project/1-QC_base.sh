#!/usr/bin/env bash
#SBATCH --mem=10GB
#SBATCH --job-name="cleanbase"
#SBATCH --output=1-QC_base.log


#########
# $1 - the base dataset (GWAS summary statistics)
# $2 - project directory, where the base data and scripts are stored
########


# ===> CHANGE HERE =================================================================================

# BASEFILE=$1
PROJECTDIR="/home/s.defina/EGG" # project folder
FILES="CAD T2D"  # "Inouye Mahajan"

# ==================================================================================================

# Set working directory to the project folder
cd $PROJECTDIR

# ==================================================================================================
echo " "
echo "======================== CLEANING BASE FILE(S) ========================"
echo " "

for base in $FILES; do

# Visualize
zcat PGS-${base}.txt.gz | head -n 18

echo " "
echo "Original file size:"
zcat PGS-${base}.txt.gz | wc -l

# Duplicate SNPs
gunzip -c PGS-${base}.txt.gz | awk '{seen[$1]++; if(seen[$1]==1){ print}}' | gzip - > PGS-${base}.nodup.gz

echo "Duplicates removed. New file size:"
zcat PGS-${base}.nodup.gz | wc -l

# Ambiguous SNPs
gunzip -c PGS-${base}.nodup.gz | awk '!( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="G" && $5=="C") || ($4=="C" && $5=="G")) {print}' | gzip > PGS-${base}.QC.gz

echo "Ambiguous SNPs removed. New file size:"
zcat PGS-${base}.QC.gz | wc -l

# Remove the first 14 rows (metadata)
zcat PGS-${base}.QC.gz | tail -n +15 > PGS-${base}.clean.txt

echo "Final number of rows:"
cat PGS-${base}.clean.txt | wc -l

# Add chr:pos column 
awk -F'\t' '{$1=++i FS $1"\t"$2":"$3":"$4":"$5;}1' OFS='\t' PGS-${base}.clean.txt | cut -f2- > PGS-${base}.QC.txt

rm PGS-${base}.nodup.gz
rm PGS-${base}.QC.gz
rm PGS-${base}.clean.txt

# No overlap file # DOES NOT WORK (??)
# grep -vf overlapped_snp_list.txt PGS-${base}.QC.txt > PGS-${base}.no_overlap.txt

# echo " "
# echo "Overlapping SNPs removed. 'no_overlap' file size:"
# cat PGS-${base}.no_overlap.txt | wc -l


echo " "
cat PGS-${base}.QC.txt | head -n 5

echo " "
echo " "

done


