
projdir = "/home/s.defina/EGG/"
setwd(projdir)


models <- expand.grid(prscore = c("CAD","T2D"),
                      parent = c("moms","dads"),
                      outcome = c("bw","pw") #,
                      # snp_ovlp <- c("yes","no")
                      # ethnic <- c("eur","noneur")
                      
)

run_regression <- function(prsscore, parent, outcome, ethn = "eur"){
    
    # Read in files
    prs_score <- read.table(paste0("./results/PRS_",prsscore,"_",ethn,".sscore"),
                            col.names = c("FID","IID","ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM","SCORE1_AVG","SCORE1_SUM"))
    
    pheno <- read.csv(paste0("./pheno/",parent,"_",ethn,"_",outcome,".csv"))
    # Merge them
    data <- merge(pheno, prs_score, by.x="PARENT", by.y="IID", all.x=TRUE)
    
    cat(paste(prsscore, "~", toupper(outcome), "in", ethn, parent, 
    "----------------------------------------\n"))
    
    cat("\nN =", nrow(data),"\n\n")
    
    if (outcome == "bw") { outc = "WEIGHT" } else { outc = "PLAWGHT"}
    
    cat(toupper(outcome), "summary:\n")
    print(summary(data[,outc]))
    cat(paste0("\n",toupper(outcome)), "SD =", round(sd(data[,outc], na.rm=TRUE), 3))
    
    cat("\n\nPRS score summary:\n")
    print(summary(data$SCORE1_AVG))
    cat("\nPRS score SD =", round(sd(data$SCORE1_AVG, na.rm=TRUE), 3), "\n")
    
    # Compute z-scores
    data$PRS_zscore <- (data$SCORE1_AVG - mean(data$SCORE1_AVG, na.rm=TRUE)) / sd(data$SCORE1_AVG, na.rm=TRUE)
    # print(summary(data$PRS_zscore))
    
    data[,paste0(outcome,"_zscore")] <- (data[,outc] - mean(data[,outc], na.rm=TRUE)) / sd(data[,outc], na.rm=TRUE)
    
    # Define regression 
    f = paste0(outcome,"_zscore ~ PRS_zscore + GENDER + GESTBIR + ",paste0("C",1:10, collapse = " + "))
    # fit it 
    mod = lm(as.formula(f), data=data)
    
    cat("\nformula: ", f, "\n")
    # save output
    print(summary(mod))
    # save(mod, file = paste0("./results/",prsscore,"/",prsscore,"_",outcome,"_",parent,"_",ethn,".rds"))
}

sink(file="./4-Regressions.log")
allmods <- do.call(mapply, c(run_regression, unname(models)))
sink()