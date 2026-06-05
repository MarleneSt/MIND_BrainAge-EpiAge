
#VB: 
# - Replaced datMeth with DNAm.
# - Fixed the typo in coef vs coefs.
# - names(PCBrainAge_prediction) <- rownames(DNAm) # VB: added this to keep sample names

#' calcPCBrainAge
#'
#' @description A function to calculate the predictor of PCBrainAge
#'
#' @param DNAm A matrix of methylation beta values. Rows = samples, columns = CpGs. Column names are required.
#' @param pheno Optional: Sample phenotype data (also with samples as rows). The calculated PCBrainAge will be appended to this.
#' @param CpGImputation An optional named vector with mean values for each CpG. If not provided, the function uses `imputeMissingBrainCpGs`.
#'
#' @return If `pheno` is provided, returns `pheno` with a new `PCBrainAge` column. Otherwise, returns a vector of predicted values.
#' @export
#'
#' @examples
#' calcPCBrainAge(exampleBetas, examplePheno)
calcPCBrainAge <- function(DNAm, pheno = NULL, CpGImputation = NULL) {
  
  ##########################
  ### Load Model Objects ###
  ##########################
  
  # # Option 1 (Preferred): Load pre-built model from package data
  # data("PCBrainAge_Model", envir = environment())
  
  # Uncomment and use below only if you are *not* loading the pre-built object
  PCBrainAge_Model <- list()
  PCBrainAge_Model$rotation <- rbind(rotation1, rotation2, rotation3)
  PCBrainAge_Model$center <- c(centering1, centering2)
  PCBrainAge_Model$coefs <- modelFit$coefs
  PCBrainAge_Model$intercept <- modelFit$intercept

  ############################################
  ### Ensure All Required CpGs Are Present ###
  ############################################
  
  requiredCpGs <- rownames(PCBrainAge_Model$rotation)
  missingCpGs <- setdiff(requiredCpGs, colnames(DNAm))
  
  if (length(missingCpGs) > 0) {
    DNAm[, missingCpGs] <- NA  # Add missing CpGs as NA columns
    
    if (is.null(CpGImputation)) {
      data("imputeMissingBrainCpGs", envir = environment())
      CpGImputation <- imputeMissingBrainCpGs
    }
    
    # Fill missing CpGs using vectorized imputation
    DNAm[, missingCpGs] <- matrix(CpGImputation[missingCpGs],
                                  nrow = nrow(DNAm),
                                  ncol = length(missingCpGs),
                                  byrow = TRUE)
    
    message("Missing CpGs were imputed using provided or default mean values.")
  }
  
  ##################################
  ### Reorder and Impute Missing ###
  ##################################
  
  DNAm <- DNAm[, requiredCpGs]  # Ensure correct order
  
  if (any(is.na(DNAm))) {
    DNAm <- apply(DNAm, 2, meanImpute)
    message("Mean imputation applied to missing values within DNAm.")
  }
  
  ##################################
  ### Final Prediction Step ###
  ##################################
  
  centered <- sweep(as.matrix(DNAm), 2, PCBrainAge_Model$center, "-")
  projection <- centered %*% PCBrainAge_Model$rotation
  PCBrainAge_prediction <- as.vector(projection %*% PCBrainAge_Model$coefs + PCBrainAge_Model$intercept)
  names(PCBrainAge_prediction) <- rownames(DNAm) # VB: added this to keep sample names
  
  
  message("PCBrainAge successfully calculated!")
  
  ############################
  ### Return Final Result ###
  ############################
  
  if (is.null(pheno)) {
    out <- data.frame(
      Sample_Name = rownames(DNAm),
      PCBrainAge = PCBrainAge_prediction
    )
  } else {
    pheno$PCBrainAge <- PCBrainAge_prediction
    out <- pheno
  }
  return(out)
}
