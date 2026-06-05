# needs these libraries installed
# library(impute)
# library(matrixStats)
# library(meffil)
# library(maxprobes) # for cross-reactive probes

# _______________________________________________#
#       POSTPROCESS DATA: 450K and epic 
# _______________________________________________#

preprocess_methylation <- function(meth_matrix,
                                   remove_sexchr = FALSE,
                                   remove_crossreactive = FALSE,
                                   array_type = c("450k", "EPIC"),
                                   output_prefix = "DNAm.imputed_winsorized",
                                   time_point = NULL,
                                   epic_chunk_size = 100000,
                                   log_file = "preprocess_log.txt") {
  
  array_type <- match.arg(array_type)
  
  # ------------------------
  # Determine base folder & file prefix
  # ------------------------
  #base_folder <- ifelse(dirname(output_prefix) == ".", getwd(), dirname(output_prefix))
  base_folder <- ifelse(dirname(output_prefix) == ".", getwd(), dirname(output_prefix))
  file_prefix <- tools::file_path_sans_ext(basename(output_prefix))
  
  # ------------------------
  # Folders
  # ------------------------
  descriptives_folder <- file.path(base_folder, "descriptives")
  logs_folder        <- file.path(base_folder, "logs")
  imputed_folder     <- file.path(base_folder, "imputed")
  
  for (f in c(descriptives_folder, logs_folder, imputed_folder)) {
    if (!dir.exists(f)) dir.create(f, recursive = TRUE)
  }
  
  # ----------------------------------------------------------
  # Capture all messages in a log file 
  # ----------------------------------------------------------
  
  con <- file(log_file, open = "a")  # append mode
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    close(con)
  }, add = TRUE)
  
  message("\n--- Preprocessing started at ", Sys.time(), " for array: ", array_type, " ---")
  
  # ----------------------------------------------------------------------
  # STEP 0a: REMOVE SEX PROBES (optional)
  # ----------------------------------------------------------------------
  
  # Step 0a: Remove sex chromosome probes if requested (using meffil)
  if (remove_sexchr) {
    message("Removing probes on sex chromosomes with meffil...")
    # get the features for your array type (e.g., "450k" or "epic")
    featureset <- tolower(array_type)  # just in case, normalise case
    features <- meffil.get.features(featureset)
    
    # subset to autosomal probes only
    autosomal_probes <- features$name[features$chr %in% paste0("chr", 1:22)]
    message("Probes before sex chromosome removal: ", nrow(meth_matrix))
    
    # subset methylation matrix to autosomal probes only
    meth_matrix <- meth_matrix[rownames(meth_matrix) %in% autosomal_probes, , drop = FALSE]
    message("Probes after sex chromosome removal: ", nrow(meth_matrix))
  }
  
  # ----------------------------------------------------------------------
  # STEP 0b: REMOVE CROSS-REACTIVE PROBES (optional)
  # ----------------------------------------------------------------------
  
  # Step 0b: Remove cross-reactive probes if requested
  if (remove_crossreactive) {
    message("Removing cross-reactive probes...")
    xloci_450K <- unique(maxprobes::xreactive_probes(array_type = "450K"))
    xloci_EPIC <- unique(unlist(maxprobes::xreactive_probes(array_type = "EPIC")))
    xloci <- unique(c(xloci_450K, xloci_EPIC))
    
    probes_before <- nrow(meth_matrix)
    meth_matrix <- meth_matrix[!(rownames(meth_matrix) %in% xloci), , drop = FALSE]
    probes_after <- nrow(meth_matrix)
    
    message("Cross-reactive probes removed: ", probes_before - probes_after)
    message("Remaining probes after cross-reactive probe removal: ", probes_after)
  }
  
  # ----------------------------------------------------------------------
  # STEP 1: Filter out probes with >30% missing
  # ----------------------------------------------------------------------
  nas.per.cpg <- rowSums(is.na(meth_matrix))
  number.samples <- ncol(meth_matrix)
  is.cpg.more.than.thirtyPercent.is.na <- (nas.per.cpg / number.samples) > 0.30
  meth_filtered <- meth_matrix[!is.cpg.more.than.thirtyPercent.is.na, , drop = FALSE]
  
  message("Filtered out probes with >30% missing values: ", nrow(meth_matrix) - nrow(meth_filtered))
  message("Remaining probes after probe filtering: ", nrow(meth_filtered), ", samples: ", ncol(meth_filtered))
  message("Missing values before imputation: ", sum(is.na(meth_filtered)))
  
  # ----------------------------------------------------------------------
  # STEP 2: Order probes by chromosome and position before imputation
  # ----------------------------------------------------------------------
  message("Ordering probes by genomic position...")
  
  featureset <- tolower(array_type)  # e.g. "450k" or "epic"
  annotation_data <- meffil.get.features(featureset)
  
  # remove probes with missing chr or pos info
  annotation_data <- annotation_data[!is.na(annotation_data$chr) & !is.na(annotation_data$pos), ]
  
  # create numeric chromosome column for ordering
  summary(as.factor((annotation_data$chromosome)))
  annotation_data$chromosome[annotation_data$chr == "chrX"] <- 23
  annotation_data$chromosome[annotation_data$chr == "chrY"] <- 24
  annotation_data$chr_num <- as.numeric(gsub("chr", "", annotation_data$chromosome))
  
  # order by chromosome then position
  annotation_data <- annotation_data[order(annotation_data$chr_num, annotation_data$pos), ]
  
  # subset to probes in methylation matrix, in this order
  ordered_probes <- annotation_data$name
  ordered_probes <- ordered_probes[ordered_probes %in% rownames(meth_filtered)]
  
  # reorder methylation matrix rows accordingly
  meth_filtered <- meth_filtered[ordered_probes, , drop = FALSE]
  message("Remaining probes after ordering by chr and pos: ", nrow(meth_filtered), ", samples: ", ncol(meth_filtered))
  
  # ----------------------------------------------------------------------
  # STEP 3: Remove samples (individuals) with >30% missingness
  # ----------------------------------------------------------------------
  sample_na_prop <- colMeans(is.na(meth_filtered))
  high_missing_samples <- names(sample_na_prop[sample_na_prop > 0.30])
  
  if (length(high_missing_samples) > 0) {
    message("Removing ", length(high_missing_samples), " samples with >30% missing values.")
    meth_filtered <- meth_filtered[, !(colnames(meth_filtered) %in% high_missing_samples), drop = FALSE]
    message("Remaining probes: ", nrow(meth_filtered), ", samples: ", ncol(meth_filtered))
  } else {
    message("No samples with >30% missing values.")
    message("Remaining probes: ", nrow(meth_filtered), ", samples: ", ncol(meth_filtered))
    
  }
  
  # ----------------------------------------------------------------------
  # STEP 4: Impute missing values using knn
  # ----------------------------------------------------------------------
  # Use chunking for EPIC if needed (due to memory)
  message("Imputing missing values using knn...")
  if (array_type == "EPIC") {
    message("Using chunked imputation for EPIC array with chunk size = ", epic_chunk_size)
    n_probes <- nrow(meth_filtered)
    chunks <- split(seq_len(n_probes), ceiling(seq_len(n_probes) / epic_chunk_size))
    
    imputed_chunks <- lapply(chunks, function(idx) {
      message("Imputing chunk with probes ", min(idx), " to ", max(idx))
      impute.knn(meth_filtered[idx, , drop = FALSE])$data
    })
    meth_imputed <- do.call(rbind, imputed_chunks)
    rownames(meth_imputed) <- rownames(meth_filtered)
    colnames(meth_imputed) <- colnames(meth_filtered)
  } else {
    imputed <- impute.knn(meth_filtered)
    meth_imputed <- imputed$data
  }
  
  message("Missing after imputation: ", sum(is.na(meth_imputed)))
  message("Value range before imputation: ", paste(range(meth_filtered, na.rm = TRUE), collapse = " - "))
  message("Value range after imputation: ", paste(range(meth_imputed), collapse = " - "))
  
  message("Mean before imputation: ", paste(mean(meth_filtered, na.rm = TRUE), collapse = " - "))
  message("Mean after imputation: ", paste(mean(meth_imputed), collapse = " - "))
  
  # ----------------------------------------------------------------------
  # STEP 5: Winsorize extreme values (>3*IQR per probe)
  # ----------------------------------------------------------------------
  message("Winsorizing extreme values...")
  
  q1 <- rowQuantiles(meth_imputed, probs = 0.25, na.rm = TRUE)
  q3 <- rowQuantiles(meth_imputed, probs = 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  
  lower <- q1 - 3 * iqr
  upper <- q3 + 3 * iqr
  
  meth_winsor <- meth_imputed
  
  # For tracking
  winsorized_cpg <- logical(nrow(meth_winsor))
  num_values_replaced <- 0
  
  for (i in seq_len(nrow(meth_winsor))) {
    original_values <- meth_winsor[i, ]
    winsorized_values <- pmin(pmax(original_values, lower[i]), upper[i])
    
    # Count how many values were changed
    values_changed <- sum(winsorized_values != original_values, na.rm = TRUE)
    num_values_replaced <- num_values_replaced + values_changed
    
    # Flag if at least one was changed
    winsorized_cpg[i] <- values_changed > 0
    
    # Store back
    meth_winsor[i, ] <- winsorized_values
  }
  
  # Final checks
  stopifnot(all(rownames(meth_winsor) == rownames(meth_imputed)))
  stopifnot(all(colnames(meth_winsor) == colnames(meth_imputed)))
  
  # Summaries
  num_cpgs_winsorized <- sum(winsorized_cpg, na.rm = TRUE)
  pct_cpgs_winsorized <- 100 * num_cpgs_winsorized / nrow(meth_winsor)
  
  total_values <- prod(dim(meth_winsor))
  pct_values_replaced <- 100 * num_values_replaced / total_values
  
  message(sprintf("Number of CpGs with at least one value winsorized: %d (%.2f%%)", 
                  num_cpgs_winsorized, pct_cpgs_winsorized))
  message(sprintf("Total values winsorized: %d (%.4f%% of all data points)", 
                  num_values_replaced, pct_values_replaced))
  
  # check range and mean after winsorization
  message("Value range after winsorization: ", paste(range(meth_winsor), collapse = " - "))
  message("Mean after winsorization: ", mean(meth_winsor))
  
  # Save names of affected CpGs
  winsorized_cpg_names <- rownames(meth_winsor)[winsorized_cpg]
  #writeLines(winsorized_cpg_names, paste0(array_type, "_", time_point, "_winsorized_cpg_names.txt"))
  #writeLines(winsorized_cpg_names, paste0(file_prefix, "_", time_point, "_winsorized_cpg_names.txt"))
  writeLines(winsorized_cpg_names, file.path(descriptives_folder, paste0(file_prefix, "_", time_point, "_", array_type, "_winsorized_cpg_names.txt"))
  )
  
  
  # ----------------------------------------------------------------------
  # STEP 6: Compute per-CpG descriptives before and after preprocessing
  # ----------------------------------------------------------------------
  #source("functions/cpg_descriptives.R")
  message("Computing per-CpG descriptives before and after preprocessing...")
  
  cpg_descriptives_before <- cpg_descriptives(meth_matrix)  # before filtering/imputation/winsor
  cpg_descriptives_after  <- cpg_descriptives(meth_winsor)   # after imputation + winsorization
  
  aligned <- merge(cpg_descriptives_before, cpg_descriptives_after, by = "CpG",
                   suffixes = c("_before", "_after"), all = TRUE)
  
  #descriptives_file <- paste0(output_prefix, "_", time_point, "_cpg_level_descriptives.rds")
  descriptives_file <- file.path(descriptives_folder, paste0(file_prefix, "_", time_point, "_", array_type, "_cpg_level_descriptives.rds"))
  saveRDS(aligned, file = descriptives_file)
  message("Per-CpG descriptives saved to: ", descriptives_file)

  # ----------------------------------------------------------------------
  # STEP 7: Save imputed + winsorized data
  # ----------------------------------------------------------------------
  saveRDS(as.matrix(meth_winsor), file.path(imputed_folder, paste0(file_prefix, "_", time_point, "_", array_type, ".rds")))
  message("Preprocessing complete. Data saved to: ", file.path(imputed_folder, paste0(file_prefix, "_", time_point, "_", array_type, ".rds")))
  
  # ----------------------------------------------------------------------
  # TIDYING: Move log file to 'logs' folder
  # ----------------------------------------------------------------------
  if (file.exists(log_file)) {
    file.rename(log_file, file.path(logs_folder, basename(log_file)))
  }
  
  message("--- Preprocessing finished at ", Sys.time(), " ---\n")
  
  return(as.matrix(meth_winsor))
}
