# ------------------------------------------------------------------------------
# 0. Prepare phenotype data ----------------------------------------------------
# ------------------------------------------------------------------------------

require(foreign)

# Point to the project folder
# if (exists("fileloc") == F) { fileloc <- file.choose() }
# datadir <- dirname(fileloc)

datadir <- '/home/s.defina/EGG/pheno/raw/'

# list.files(datadir)

readit <- function(filename, sel=c()){
  out <- suppressWarnings(foreign::read.spss(file.path(datadir,paste0(filename,'.sav')), to.data.frame=TRUE))
  if (length(sel)>0){ out <- out[,sel] }
  # print(summary(out))
  return(out)
}

sink(file.path(datadir, '../../0-pheno_prep.log'))

# Load GWA sample
geno <- readit('SelectionGenR-Parents_PARENT_28022023',c('PARENT','Type'))
nonE <- readit('GENR_Parents_NonEurSamples_28022023',c('PARENT','NonEuropean'))
# Merge
geno <- merge(geno, nonE, by='PARENT', all=TRUE)
# Replace NA with 0s
geno$NonEuropean <- as.factor(ifelse(is.na(geno$NonEuropean), 0, 1))

# Add genetic PCs
PCA <- readit('GenRParents_Indep_marker_PCs_18072023', c('Parent',paste0('C',1:10)))
geno <- merge(geno, PCA, by.x='PARENT',by.y='Parent', all=TRUE)

cat("\nGENETIC DATA FILE\n")
summary(geno) # Note: 7236 mothers and 4440 partners
# any(duplicated(geno$PARENT)) # FALSE

# Load phenotype and exclusion criteria
phen1 <- readit('CHILD-ALLGENERALDATA_24102022', c('IDM','IDC','MOTHER','PARTNER','PARTBF','PARITY','TWIN','GENDER','GESTBIR','WEIGHT'))
phen2 <- readit('FETALPLACENTA_22112016', c('IDM','PLAWGHT'))
# Merge
phen <- merge(phen1, phen2, by='IDM', all.x=TRUE)

cat("\nPHENOTYPE DATA FILE\n")
summary(phen)

# Perform sample selection
clean_samp <- function(parent, eur=0) { 
  cat('\n=============', parent, '=============\n')
  dset <- merge(geno, phen, by.x='PARENT', by.y=parent, all=FALSE)
  n = nrow(dset)
  
  
  cat(sum(duplicated(dset$PARENT)), 'duplicated parent IDs. (Orig:', 
      n, '-', n-sum(duplicated(dset$PARENT)), ')')
  
  cat('\n\nOnly term births (37-42 weeks): \n')
  n = nrow(dset)
  dset <- dset[dset$GESTBIR > 37 & dset$GESTBIR < 42,]
  cat(n, '-', n-nrow(dset), '=', nrow(dset) )
  
  cat('\n\nOnly non-twin first-born: \n')
  n = nrow(dset)
  dset <- dset[(!is.na(dset$PARITY) & dset$PARITY==0) & 
               (!is.na(dset$TWIN) & dset$TWIN=='No'),]
  cat(n, '-', n-nrow(dset), '=', nrow(dset) )
  
  cat('\n\nOnly (non-)/Europeans: \n') # TODO: supplementary analysis non-euro
  n = nrow(dset)
  dset <- dset[dset$NonEuropean==eur,]
  cat(n, '-', n-nrow(dset), '=', nrow(dset) )
  
  cat('\n\nAvailable birthweight and placental weight: \n')
  dsetW <- dset[!is.na(dset$WEIGHT),]
  cat('B.weight: ',n, '-', n-nrow(dsetW), '=', nrow(dsetW),'\n' )
  dsetP <- dset[!is.na(dset$PLAWGHT),]
  cat('Placent.: ', n, '-', n-nrow(dsetP), '=', nrow(dsetP))
  
  # Correct issue with duplicate parent 
  dsetW <- dsetW[-which(duplicated(dsetW$PARENT)),]
  
  id_list <- merge(dsetW, dsetP, by='PARENT', all=TRUE)['PARENT']
  
  cat('\n\nRemaining duplicates: ', sum(duplicated(dsetW$PARENT)), sum(duplicated(dsetP$PARENT)), '\n')
  
  # Save files to csv
  if (parent=='MOTHER') { name = 'moms' } else { name = 'dads'}
  if (eur==0) { ethn = 'eur' } else { ethn = 'noneur'}
  write.csv(dsetW, file.path(datadir, paste0('../',name,'_',ethn,'_bw.csv')), row.names = FALSE)
  write.csv(dsetP, file.path(datadir, paste0('../',name,'_',ethn,'_pw.csv')), row.names = FALSE)
  
  return(id_list) # list(dsetW,dsetP,id_list)
}

cat('\n==================================\nEUROPEANS\n==================================\n')
moms = clean_samp('MOTHER')
dads = clean_samp('PARTNER')
cat('\n==================================\nNON-EUROPEANS\n==================================\n')
moms_noneur = clean_samp('MOTHER', eur=1)
dads_noneur = clean_samp('PARTNER', eur=1)

sink()

# Save participant list 
save_plink_format <- function(x, filename){
  cat(nrow(x), ' IDs.\n')
  plink_form <- data.frame(x, x)
  names(plink_form) <- c('FID','IID')
  write.table(x = plink_form, 
              file = file.path(dirname(datadir), filename),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
}


save_plink_format(rbind(moms,dads), 'id_list_eur.txt')
save_plink_format(rbind(moms_noneur,dads_noneur), 'id_list_noneur.txt')
