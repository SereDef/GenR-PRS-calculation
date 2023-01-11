library(data.table)
library(stringr)

data.table::setDTthreads(10) # set number of threds for parallel functions 
args = commandArgs(trailingOnly = TRUE) # provides access to a copy of the command line arguments supplied

if (length(args) == 0) {
  stop("Supply DOSAGES file as well as SNPs file")
} else {
  DOS <- args[1]
  SNP <- args[2]
}

cat("=====================================\nReading in DOSAGE and SNP/GWAS files")
DOSAGES <- data.table::fread(file = DOS, header = TRUE, stringsAsFactors = FALSE)
DOSAGES[, SNP := paste(CHR, POS, sep = ":")]

# Separate AF, MAF and imputation quality in separate columns 
info <- stringr::str_split_fixed(DOSAGES$INFO, ';', 4)

DOSAGES$AF  <- as.numeric(stringr::str_split_fixed(info[,1], '=', 2)[,2])
DOSAGES$MAF <- as.numeric(stringr::str_split_fixed(info[,2], '=', 2)[,2])
DOSAGES$R2  <- as.numeric(stringr::str_split_fixed(info[,3], '=', 2)[,2])

DOSAGES <- DOSAGES[,-6] # which(names(DOSAGES)=="INFO")

SNPS <- data.table::fread(file = SNP, header = TRUE, stringsAsFactors = FALSE)

cat("\n-------------------------------------\nHarmonizing files\n")
# order by chr:position
SNPS <- SNPS[order(SNPS$SNP)]; DOSAGES <- DOSAGES[order(DOSAGES$SNP)]

# There should not be any variant in the target file that is not in the base file 
if(!identical(which(!DOSAGES$SNP %in% SNPS$SNP),integer(0))) { 
  stop("Target file contains variants that are not in the base data. Check the selection process.\n") }

# Identify SNPs missing in the target data
MIS <- SNPS[which(!SNPS$SNP %in% DOSAGES$SNP)]$SNP
if (!identical(MIS, character(0))) { 
  cat(length(c(MIS)), "SNP(s) missing in target data: ", MIS, "\nRemoving them.")
  # Remove un-matched SNPs 
  SNPS <- SNPS[which(SNPS$SNP %in% DOSAGES$SNP)]
} else { cat("No missing SNPs in target data (yey).") }

cat("\n-------------------------------------\nRemoving duplicated/multiallelic SNPs\n")
# Identify multiallelic SNPs
MUL <- DOSAGES[!DOSAGES$ALT %in% c('A','C','G','T')]$SNP
if (!identical(MUL, character(0))) { 
  cat(length(c(MUL)), "multiallelic SNP(s) in target data: ", MUL, "\nRemoving them.\n")
  DOSAGES <- DOSAGES[!DOSAGES$SNP %in% c(MUL),]
} else { cat("No multiallelic SNPs detected in target data.\n") }

table(DOSAGES[,c("ALT","REF")])

# Identify ambiguous SNPs 
AMB <- DOSAGES[(DOSAGES$ALT=="A" & DOSAGES$REF=="T")|(DOSAGES$ALT=="T" & DOSAGES$REF=="A")| 
               (DOSAGES$ALT=="C" & DOSAGES$REF=="G")|(DOSAGES$ALT=="G" & DOSAGES$REF=="C"), c('SNP','ALT','REF')]
if (nrow(AMB)>0) {
  at = AMB[AMB$ALT=='A']$SNP; ta = AMB[AMB$ALT=='T']$SNP; cg = AMB[AMB$ALT=='C']$SNP; gc = AMB[AMB$ALT=='G']$SNP
  cat("\n", nrow(AMB), "ambiguous SNPs in target data:", length(at),"A-T pairs [",at,"];",
      length(ta),"T-A pairs [",ta,"];",length(cg),"C-G pairs [",cg,"]; and",length(gc),"C-G pairs [",gc,"].",
      "\nRemoving them.\n")
  DOSAGES <- DOSAGES[!DOSAGES$SNP %in% c(AMB$SNP),]
  # table(AMB[,c("ALT","REF")])
}

# Identify duplicated SNPs
DUP <- DOSAGES[which(duplicated(DOSAGES$SNP))]$SNP
if (!identical(DUP, character(0))) { 
  cat("\n", length(c(DUP)), "duplicated SNP(s) in target data: ", DUP, "\nRemoving them.\n")
  DOSAGES <- DOSAGES[!DOSAGES$SNP %in% c(DUP),] # <==== this removes everything (change if not desired)
} else { cat("No duplicated SNPs detected in target data.\n") }

# Update the base dataset
SNPS <- SNPS[which(SNPS$SNP %in% DOSAGES$SNP)]

cat("\n-------------------------------------\nWorking with", dim(DOSAGES)[1], "SNPs.")

if (!identical(SNPS$SNP, DOSAGES$SNP)) {
  stop("\nERROR: target and base data are not aligned!")
}

cat("\n-------------------------------------\nInverting beta estimates\n")
count <- 0; DROPPED <- c()
for (i in 1:nrow(SNPS)) {
  if (SNPS[i, "EA"] != DOSAGES[i, "ALT"]) { count <- count + 1
    if (SNPS[i, "EA"] == DOSAGES[i, "REF"]) {
      # Invert beta estimate in the base data
      SNPS[i, 'beta'] <- -SNPS[i, 'beta']
    } else { cat("Problem with phasing in SNP(s)", DOSAGES[i, SNP], "; Dropping them.\n")
      DROPPED <- c(DROPPED, i)
    }
  }
}
cat(count, "beta estimates were inverted in the base dataset to match target data.")

# cat("\n-------------------------------------\nMatching target data with base effect alleles\n")
# count <- 0; DROPPED <- c()
# for (i in 1:nrow(SNPS)) {
#   if (SNPS[i, "EA"] != DOSAGES[i, "ALT"]) { count <- count + 1
#     if (SNPS[i, "EA"] == DOSAGES[i, "REF"]) {
#       # Swap effect and reference alleles in target data
#       DOSAGES[i, "REF"] <- DOSAGES[i, "ALT"]
#       DOSAGES[i, "ALT"] <- SNPS[i, "EA"]
#       # correct individual dosage data to allign with new effect/reference
#       DOSAGES[i, 6:ncol(DOSAGES)] <- -1 * (DOSAGES[i, 6:ncol(DOSAGES)] - 2)
#     } else { cat("Problem with phasing in SNP(s)", DOSAGES[i, SNP], "; Dropping them.\n")
#       DROPPED <- c(DROPPED, i)
#     }
#   }
# }
# cat(count, "allele pairs were inverted in the target dataset to match base data.")

# if (!identical(DOSAGES$ALT, SNPS$EA)) { stop("\nERROR: target and base data not harmonized!") } 

# cat("\n-------------------------------------\nChanging beta signs to align with effect increasing alleles.\n")
# count <- 0
# for (i in 1:nrow(SNPS)) {
#   if (SNPS[i, 'beta'] < 0) { count <- count + 1
#     SNPS[i, beta := -1 * beta] # invert effect estimate 
#     SNPS[i, EAF := 1 - EAF] # invert effect allele frequency
#     # Swap effect and reference alleles in target and base data
#     DOSAGES[i, ALT := REF] 
#     DOSAGES[i, ]$REF <- SNPS[i, ]$EA
#     SNPS[i, ]$EA <- DOSAGES[i, ]$ALT
#     # correct individual dosage data to allign with new effect/reference
#     DOSAGES[i, 6:ncol(DOSAGES)] <- -1 * (DOSAGES[i, 6:ncol(DOSAGES)] - 2)
#   }
# }
# cat("All beta estimates are positive. (", count, " estimates were inverted).\n", sep="") 

# if (any(SNPS$BETA < 0)) { stop("\nERROR: some effect sizes still negative.\n") }

# if (!identical(DOSAGES$ALT, SNPS$EA)) { stop("\nERROR: target and base data not harmonized!")
# } else { cat("\n-------------------------------------\nAll alleles seem harmonized.") }
cat("\n-------------------------------------\nFiltering based on MAF and imputation quality (R2).\n")

rare <- DOSAGES[DOSAGES$MAF < 0.01, ]
if (nrow(rare) > 0) {
  cat(nrow(rare), 'variants with MAF < 0.01: ', rare$SNP,'\nRemoving them.\n')
  DOSAGES <- DOSAGES[DOSAGES$MAF >= 0.01, ]
} else {
  cat('All variants have MAF > 0.01.\n')
}

impq <- DOSAGES[DOSAGES$R2 < 0.5, ]
if (nrow(impq) > 0) {
  cat(nrow(impq), 'variants with imputation quality (R2) < 0.5: ', impq$SNP,'\nRemoving them.\n')
  DOSAGES <- DOSAGES[DOSAGES$R2 >= 0.5, ]
} else {
  cat('All variants have imputation quality (R2) > 0.5.\n')
}

# Update the base dataset
SNPS <- SNPS[which(SNPS$SNP %in% DOSAGES$SNP)]

# Save list of excluded SNPs
excl <- c(MIS, MUL, AMB$SNP, DUP, rare$SNP, impq$SNP)
write.table(x = excl, file = file.path(dirname(SNP), "EXCLUDED-SNPS-GENR.txt"), 
            quote = F, row.names = F, col.names = F, sep = "\n") 

cat("\n-------------------------------------\nExtracting PGS based on", dim(DOSAGES)[1], "SNPs.\n")

PRS <- DOSAGES[, 6:ncol(DOSAGES)]
WPRS <- as.data.frame(as.matrix(PRS) * SNPS$beta)

final_df <- data.frame(IDC = names(PRS), PGS = colSums(PRS), wPGS = colSums(WPRS))
final_df$wPGS_scaled <- scale(final_df$wPGS)

write.table(x = final_df, file = file.path(dirname(SNP), "PGS-GENR.txt"), 
            quote = F, row.names = F, col.names = T, sep = "\t")
cat("Done.\n-------------------------------------\n")
