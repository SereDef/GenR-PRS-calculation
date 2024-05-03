
args = commandArgs(trailingOnly = TRUE) # provides access to a copy of the command line arguments supplied

if (length(args) == 0) { stop("Supply het file location!") } else { path <- args[1] }

setwd(path)

exclude <- character() # initialize

for (i in c(1:22)) {
  het <- read.table(dir(pattern = "chr.*het")[i], header = T)
  het$HET_RATE <- (het$N_SITES - het$O.HOM.)/het$N_SITES
  avg <- mean(het$HET_RATE)
  sdev <- sd(het$HET_RATE)
  exclude <- c(exclude, het[which(het$HET_RATE > 3*sdev + avg | het$HET_RATE < avg - 3*sdev), "INDV"])
}

exclude <- unique(exclude)
extable <- data.frame("PIF" = exclude, "IID" = exclude)

write(exclude, "excl_het.txt")
write.table(extable, "excl_het_tab.txt", quote = F, row.names = F, sep = " ")

#het <- read.table("output.het", head=TRUE)
#pdf("heterozygosity.pdf")
#het$HET_RATE = (het$N_SITES - het$O.HOM.)/het$N_SITES
#avg <- mean(het$HET_RATE)
#sdev <- sd(het$HET_RATE)
#exclude <- het[which(het$HET_RATE > 3*sdev + avg | het$HET_RATE < avg - 3*sdev), "INDV"]
#write(exclude, "excl_het.txt")
#hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
#dev.off()
