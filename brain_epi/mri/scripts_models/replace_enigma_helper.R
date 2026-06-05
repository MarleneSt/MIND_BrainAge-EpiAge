#===============================================================================
# replace_from_enigma(): Join target data to ENIGMA by ID and overwrite values
# - Preserves original target columns & order
# - Supports different ID column names (e.g., "ID" vs "SUBJID")
# - Modes:
#   * "overwrite_all"          : ENIGMA fully replaces target (even with NAs)
#   * "overwrite_nonmissing"   : replace only where ENIGMA is non-NA
#   * "overwrite_missing"      : fill only where TARGET is NA
#   * "overwrite_if_enigma_na" : set TARGET to NA wherever ENIGMA is NA
#
# Marlene Staginnus
#===============================================================================
replace_from_enigma <- function(
    target_df,
    enigma_df,
    id_col_target = "SUBJID",
    id_col_enigma = "SUBJID",
    mapping,                         # data.frame/tibble with cols: target, enigma
    mode = c("overwrite_if_enigma_na","overwrite_all", "overwrite_nonmissing", "overwrite_missing"),
    output_path = NULL               # optional CSV path to save
) {
  mode <- match.arg(mode)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  
  # keep original schema
  orig_names <- names(target_df)
  
  # --- key checks
  if (!id_col_target %in% names(target_df)) stop("ID column '", id_col_target, "' not found in target_df.")
  if (!id_col_enigma %in% names(enigma_df)) stop("ID column '", id_col_enigma, "' not found in enigma_df.")
  if (!all(c("target", "enigma") %in% names(mapping))) {
    stop("`mapping` must have columns named 'target' and 'enigma'.")
  }
  
  # --- work on copies; align IDs to SUBJID for join
  td <- target_df
  ed <- enigma_df
  if (id_col_target != "SUBJID") names(td)[names(td) == id_col_target] <- "SUBJID"
  if (id_col_enigma != "SUBJID") names(ed)[names(ed) == id_col_enigma] <- "SUBJID"
  
  # --- SUBJID alignment check (stop if not all target IDs are in ENIGMA)
  if (!all(td$SUBJID %in% ed$SUBJID)) {
    missing_ids <- setdiff(td$SUBJID, ed$SUBJID)
    stop("❌ Not all SUBJIDs in target are found in ENIGMA.\n",
         "Example missing IDs: ", paste(head(missing_ids, 10), collapse = ", "),
         "\n(Showing up to 10.)")
  } else {
    message("✅ All SUBJIDs in target are found in ENIGMA.")
  }
  
  # keep only mappings that actually exist in ENIGMA
  present <- mapping[mapping$enigma %in% names(ed), , drop = FALSE]
  if (nrow(present) == 0) {
    warning("No mapping columns found in ENIGMA. Returning target unchanged.")
    final <- td
  } else {
    # --- join ENIGMA subset
    joined <- dplyr::left_join(
      td,
      dplyr::select(ed, SUBJID, dplyr::all_of(present$enigma)),
      by = "SUBJID"
    )
    
    # --- overwrite logic
    for (i in seq_len(nrow(present))) {
      tgt <- present$target[i]
      src <- present$enigma[i]
      if (!src %in% names(joined) || !tgt %in% names(joined)) next
      
      # light type guard
      if (is.numeric(joined[[tgt]]) && !is.numeric(joined[[src]])) {
        joined[[src]] <- suppressWarnings(as.numeric(joined[[src]]))
      }
      
      if (mode == "overwrite_all") {
        # ENIGMA wins completely (even NA)
        joined[[tgt]] <- joined[[src]]
        
      } else if (mode == "overwrite_nonmissing") {
        # only where ENIGMA has data
        joined[[tgt]] <- dplyr::coalesce(joined[[src]], joined[[tgt]])
        
      } else if (mode == "overwrite_missing") {
        # only where TARGET is NA
        joined[[tgt]][is.na(joined[[tgt]])] <- joined[[src]][is.na(joined[[tgt]])]
        
      } else if (mode == "overwrite_if_enigma_na") {
        # set TARGET to NA wherever ENIGMA is NA
        joined[[tgt]][is.na(joined[[src]])] <- NA
      }
    }
    
    # drop ENIGMA columns we brought in
    joined <- dplyr::select(joined, -dplyr::any_of(unique(present$enigma)))
    final <- joined
  }
  
  # restore ID name if target used a different one
  if (id_col_target != "SUBJID" && "SUBJID" %in% names(final)) {
    names(final)[names(final) == "SUBJID"] <- id_col_target
  }
  
  # restore original schema/order
  final <- dplyr::select(final, dplyr::all_of(orig_names))
  
  # informative column check
  if (!identical(names(final), orig_names)) {
    missing_in_final <- setdiff(orig_names, names(final))
    extra_in_final   <- setdiff(names(final), orig_names)
    msg <- c()
    if (length(missing_in_final)) msg <- c(msg, paste0("Missing in final: ", paste(missing_in_final, collapse = ", ")))
    if (length(extra_in_final))   msg <- c(msg, paste0("Extra in final: ", paste(extra_in_final, collapse = ", ")))
    stop("❌ Column mismatch between final and original target_df.\n", paste(msg, collapse = "\n"))
  } else {
    message("✅ Column check passed: final dataset has the same columns (and order) as original target_df.")
  }
  
  if (!is.null(output_path)) {
    utils::write.csv(final, output_path, row.names = FALSE)
    message("💾 Saved: ", output_path)
  }
  
  return(final)
}
