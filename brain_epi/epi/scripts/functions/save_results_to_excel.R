# function to save output that also includes missing cpg names (long, as each model is a list, so need to reformat)
save_results_to_excel <- function(results, output_name) {
  # Convert unlisted CpGs to data frames
  cpgs_missing_df <- data.frame(
    Source = names(unlist(results$cpgs_missing)),
    CpG = unlist(results$cpgs_missing),
    row.names = NULL
  )
  
  cpgs_missing_GA_df <- data.frame(
    Source = names(unlist(results$cpgs_missing_GA)),
    CpG = unlist(results$cpgs_missing_GA),
    row.names = NULL
  )
  
  # combine output into list
  if (!is.null(results$age_GA_pred)) {
    output_list <- list(
      age_pred = results$age_pred,
      age_GA_pred = results$age_GA_pred,
      cpgs_missing = cpgs_missing_df,
      cpgs_missing_GA = cpgs_missing_GA_df
    )
  } else {
    output_list <- list(
      age_pred = results$age_pred,
      cpgs_missing = cpgs_missing_df,
      cpgs_missing_GA = cpgs_missing_GA_df
    )
  }
  
  # Write to Excel file
  writexl::write_xlsx(output_list, path = output_name)
  message("Saved results to ", output_name)
  
}
