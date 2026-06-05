# File: cpg_descriptives.R
#' Compute descriptive statistics for methylation data
#'
#' @param meth A matrix or data.frame with CpGs as rows and samples as columns.
#' @return A data.frame with CpG IDs and descriptive statistics (mean, median, sd, min, max, missing values).
#' @examples
#' cpg_descriptives(meth_cord_450k)
cpg_descriptives <- function(meth) {
  
  # ensure matrix format once
  m <- as.matrix(meth)
  
  out <- data.frame(
    CpG     = rownames(m),
    mean    = rowMeans(m, na.rm = TRUE),
    median  = matrixStats::rowMedians(m, na.rm = TRUE),
    sd      = matrixStats::rowSds(m, na.rm = TRUE),
    min     = matrixStats::rowMins(m, na.rm = TRUE),
    max     = matrixStats::rowMaxs(m, na.rm = TRUE),
    n_na    = rowSums(is.na(m)),
    prop_na = rowMeans(is.na(m)),
    stringsAsFactors = FALSE
  )
  
  return(out)
}

