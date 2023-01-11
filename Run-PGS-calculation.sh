#!/usr/bin/env bash
#SBATCH --mem=10GB

#########
# $1 - the base dataset (GWAS summary statistics)
# $2 - project directory, where the base data and scripts are stored
########


# ===> CHANGE HERE =================================================================================
BASEFILE=$1
PROJECTDIR=$2
GENDIR="/home/057600/GENR3/Imputed/1000G_PhaseIIIv5" # location of target genetic data (e.g., GENR3)

# ==================================================================================================

# Set working directory to the project folder
cd $PROJECTDIR

# Create subfolder to store results
OUTDIR="$PROJECTDIR/results"
mkdir -p $OUTDIR

# Clean the base data: rename and reorder columns and remove duplicated/multiallelic SNPs
# TODO: make interactive job work with SLURM. For now, fix the R script manually (colnames)
Rscript Clean-base.R ./$BASEFILE

INPUT=`awk -F'.' '{print $1"_QC."$2}' <<< $BASEFILE` # this is the name of the QC'd base file

echo "====================================="
awk '{if ($1 == "chr") print "Extracting SNPs from .vcf files, based on: "$1 ":" $2}' $INPUT
echo "====================================="

for i in {1..22}; do
# extract chr:position from the summary statistics file and divide them into chromosomes  
awk -v c=$i '{if ($1 == c) print $1":"$2"-"$2}' $INPUT > $OUTDIR/extract_chr${i}.txt
cat $OUTDIR/extract_chr${i}.txt | tr "\n" " " | sed 's/ $/\n/g' > $OUTDIR/snps # replace new lines with spaces 
mv $OUTDIR/snps $OUTDIR/extract_chr${i}.txt
if [ -s $OUTDIR/extract_chr${i}.txt ] # if there are variants on the chromosome (i.e. file is not empty)
then 
    echo "Chromosome: $i"
    # perform the selection
    tabix -h $GENDIR/chr${i}.dose.vcf.gz $(cat $OUTDIR/extract_chr${i}.txt) > $OUTDIR/chr${i}.recode.vcf
    # Fill in the information (person IDs and dosages)
    bcftools query -l $OUTDIR/chr${i}.recode.vcf > $OUTDIR/chr${i}.personid
    printf "SNP"$'\t'"CHR"$'\t'"POS"$'\t'"REF"$'\t'"ALT"$'\t' > $OUTDIR/chr${i}.dosages.txt
    tr "\n" "\t" < $OUTDIR/chr${i}.personid >> $OUTDIR/chr${i}.dosages.txt
    printf "\n" >> $OUTDIR/chr${i}.dosages.txt
    bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' $OUTDIR/chr${i}.recode.vcf >> $OUTDIR/chr${i}.dosages.txt
else
    echo "NOTE: No variants on chr: $i."
    rm $OUTDIR/extract_chr${i}.txt
fi
done

echo "Extraction done."
echo "-------------------------------------"

echo "Merging extracted dosages into a single file: DOSAGES-GENR.txt"
head -n 1 <$(ls -t $OUTDIR/chr*.dosages.txt | head -n 1) > ./DOSAGES-GENR.txt
tail -n +2 -q $OUTDIR/chr*.dosages.txt >> ./DOSAGES-GENR.txt
sed -i 's/\([0-9]*\)\t$/\1/g' ./DOSAGES-GENR.txt

MIS=$(comm -13 <(cut -f2,3 ./DOSAGES-GENR.txt | tail -n +2 | tr "\t" ":" | sort) <(cut -f1,2 -d";" $INPUT | tail -n +2 | tr ";" ":" | sort))

# echo "Missing SNPs: \n${MIS}"
echo $MIS > ./MISSING-SNPs-GENR.txt
echo "====================================="

echo "Cleaning and harmonizing target data"

Rscript Create-PGS.R ./DOSAGES-GENR.txt $INPUT
