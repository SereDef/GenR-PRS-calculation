library(data.table)
library(dplyr)

data.table::setDTthreads(10) # set number of threds for parallel functions 
args = commandArgs(trailingOnly = TRUE) # provides access to a copy of the command line arguments supplied

if (length(args) == 0) {
  stop("Supply SNPs file!")
} else {
  SNP <- args[1]
}

cat("\n=====================================\nCleaning of the base file")
SNPS <- data.table::fread(file = SNP, header = TRUE, stringsAsFactors = FALSE)

cat("\n=====================================\nRenaming columns\n-------------------------------------\n")
cat("Input columns:", names(SNPS), '\n')

if (!'chr' %in% names(SNPS)) {
  other_names <- c('CHR','Chr','chr_name','Chromosome')
  chr_name <- intersect(other_names, names(SNPS))
  if (identical(chr_name, character(0))) { chr_name <- readline("Indicate chromosome name column: ") }
  names(SNPS)[which(names(SNPS)==chr_name)] <- 'chr' }

if (!'pos' %in% names(SNPS)) {
  other_names <- c('POS','Pos','chr_position','Position','position')
  pos_name <- intersect(other_names, names(SNPS))
  if (identical(pos_name, character(0))) { pos_name <- readline("Indicate genomic position column: ") }
  names(SNPS)[which(names(SNPS)==pos_name)] <- 'pos' }

if (!'EA' %in% names(SNPS)) {
  other_names <- c('ea','effect_allele','Effect allele')
  ea_name <- intersect(other_names, names(SNPS))
  if (identical(ea_name, character(0))) { ea_name <- readline("Indicate effect allele column: ") }
  names(SNPS)[which(names(SNPS)==ea_name)] <- 'EA' }

if (!'beta' %in% names(SNPS)) {
  other_names <- c('BETA','effect_weight','Effect','effect')
  beta_name <- intersect(other_names, names(SNPS))
  if (identical(beta_name, character(0))) { beta_name <- readline("Indicate effect estimate column: ") }
  names(SNPS)[which(names(SNPS)==beta_name)] <- 'beta' }

# Optional: Effect allele frequency
if (!'EAF' %in% names(SNPS)) {
  other_names <- c('EA_frq','Effect allele frequency','EA_frequency','effect_allele_frequency')
  EAF_name <- intersect(other_names, names(SNPS))
  if (identical(EAF_name, character(0))) {
    EAF_present <- readline("Is there an effect allele frequency column? [y/n] ")
    if (EAF_present %in% c('y','Y','yes')) { 
      EAF_name <- readline("Indicate effect allele frequency column: ")
      names(SNPS)[which(names(SNPS)==EAF_name)] <- 'EAF' 
    } else { cat('Note: no effect allele frequency information in base data.') }
  } else { names(SNPS)[which(names(SNPS)==EAF_name)] <- 'EAF' }
}

SNPS[, SNP := paste(chr, pos, sep = ":")] # add a chr:position column

SNPS <- SNPS %>% relocate(chr, pos, SNP) # move these to the first columns (used for variant selection)

cat("New columns:", names(SNPS))
cat("\n-------------------------------------\n")
cat("NOTE: Heritability check not performed.\n")
cat("NOTE: Genome build has to correspond to GRCh37/hg19 (not checked).")
cat("\n-------------------------------------\nRemoving duplicated/multiallelic SNPs\n-------------------------------------\n")
# Order by chr:position
SNPS <- SNPS[order(SNPS$SNP)]

# Identify multiallelic SNPs
MUL <- SNPS[!SNPS$EA %in% c('A','C','G','T')]$SNP
if (!identical(MUL, character(0))) { 
  cat(length(c(MUL)), "multiallelic SNP(s) in base data: ", MUL, "\nRemoving them.")
  SNPS <- SNPS[!SNPS$SNP %in% c(MUL),]
} else { cat("No multiallelic SNPs detected.") }

table(SNPS$EA)

# Identify duplicated SNPs
DUP <- SNPS[which(duplicated(SNPS$SNP))]$SNP
if (!identical(DUP, character(0))) { 
  cat(length(c(DUP)), "duplicated SNP(s) in base data: ", DUP, "\nRemoving them.")
  SNPS <- SNPS[!SNPS$SNP %in% c(DUP),] # <==== this removes everything (change if not desired)
} else { cat("No duplicated SNPs detected.") }

# cat("\n=====================================\nRemving SNPs in high LD\n")
# ======> specific to CADSET plan! 
# to_remove <- c("rs62375246", "rs156394")

cat("\n-------------------------------------\nNOTE: LD structure not checked.\n")
if (any(SNPS$beta < 0)) {
  cat("NOTE: ", nrow(SNPS[SNPS$beta < 0,])," estimates are negative (refer to effect decreasing alleles). These will be inverted later.\n") 
} else { cat("All effect sizes are positive.") }

cat("\n-------------------------------------\nWorking with", dim(SNPS)[1], "SNPs.")

output_file_name <- gsub('\\.','_QC.',basename(SNP))

cat("\n-------------------------------------\nSaving QC'd base file: ", output_file_name, "\n")
write.table(x = SNPS, file = file.path(dirname(SNP), output_file_name), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
