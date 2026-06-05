# function to calculate epi ages using methylclock
processDNAmAge <- function(DNAm_imputed_rds, timepoint, output_dir) { 
  # Load DNAm imputed data
  MethylationData <- readRDS(file.path(output_dir, DNAm_imputed_rds))
  
  # convert to data frame if needed
  if (!is.matrix(MethylationData)) {
    message("Converting data frame to matrix")
    MethylationData <- as.matrix(MethylationData)
  }
  
  cat("Methylation data dimensions:", dim(MethylationData), "\n")
  
  # Check missing CpGs
  cpgs.missing.GA <- checkClocksGA(MethylationData)
  missCpGs <- checkClocks(MethylationData)
  
  # Estimate epigenetic age
  age.pred <- DNAmAge(MethylationData, clocks=c("Horvath", "Hannum", "skinHorvath","Levine", "PedBE",
                                                "Wu", "TL", "BLUP", "EN")) 
  print(head(age.pred))
  
  # Estimate gestational epigenetic age if birth timepoint
  if(timepoint == "birth") {
    age.pred.GA <- DNAmGA(MethylationData, clocks=c("Knight", "Bohlin", "EPIC"))
    print(head(age.pred.GA))
  }
  
  cat("Processing completed for:", DNAm_imputed_rds, "\n\n")
  
  # Return output list
  if(timepoint == "birth") {
    return(list(
      cpgs_missing_GA = cpgs.missing.GA,
      cpgs_missing = missCpGs,
      age_pred = age.pred,
      age_GA_pred = age.pred.GA
    ))
  } else {
    return(list(
      cpgs_missing_GA = cpgs.missing.GA,
      cpgs_missing = missCpGs,
      age_pred = age.pred
    ))
  }
}