# GenR-PRS-calculation
Set of scripts to calculate polygenic risk scores based on a given list of SNPs.

This set of scripts constructs a polygenic risk scores (PRS) of a certain phenotype, based on summary statistics: i.e. list of genetic variants associated with the phenotype of interest, according to an (completed) genome-wide association study (GWAS).

The script assumes you will have a list of SNPs (`SNPlist_phenotype.txt`) containing the most significant variants associated with the phenotype and their effect sizes. These summary statistics will be referred to as the **base dataset**. These files could be downloaded from the [PGS Catalogue](https://www.pgscatalog.org/) or provided by the authors of the GWAS, and they can come in many formats. You may need to accomodate the scrips a little bit to accomodate that. 

The **target dataset** will be the genetic data (i.e., `.vcf` files) available in Generation R. The sample you are typically going to work with are GenR3, GenR4 and GenRParents. 

To run the PGS calculation, simply type this command: 
```
# only if this is the first time you run this, first need to make it executable 
chmod u+x Run-PGS-calculation.sh 
# Run it (send the job to slurm)
sbatch Run-PGS-calculation.sh SNPlist_phenotype.txt /home/microsetion/progect_name
```

The scripts will perform some quality control steps on the base dataset, extract the relevant SNPs from the genetic files, perform additional QC and compute the PRSs. Should not take longer than a few minutes (depending on how many SNPs are included in your list) and it will return two files:
* `slurm-somenumber.out` : log file providing details on the QC steps
* `PGS-GENR.txt` : the unweighted and weighted PGS score for each participant (that you can then enter in your model)





