#!/usr/bin/env bash
#SBATCH --mem=10GB
#SBATCH --job-name="PRS-EGG"
#SBATCH --output=3-PRS_plink.log

# ===> CHANGE HERE =================================================================================

PROJECTDIR="/home/s.defina/EGG" # project folder (where base files are stored)
FILES="CAD T2D"  # "Inouye Mahajan"
ETHN="eur noneur"
# ==================================================================================================

# Set working directory to the project folder
cd $PROJECTDIR

# Results subfolder
OUTDIR="$PROJECTDIR/results"

# ==================================================================================================

# for ethn in $ETHN; do

for base in $FILES; do

mkdir -p $OUTDIR/${base}

plink2 \
    --bfile $OUTDIR/all_autosomal_eur \
    --score $PROJECTDIR/PGS-${base}.QC.txt 2 5 7 header ignore-dup-ids cols=+scoresums list-variants \
    --out $OUTDIR/PRS_${base}_eur
    # --q-score range_list SNP.pvaue \
    # --extract EUR.valid.snp \

done
# done

# TODO: Remove ovelapping SNPs 
# Do noneur