# --------- helpers ---------
descriptives <- function(meth, ...) {
  m <- as.matrix(meth)
  
  q1  <- matrixStats::rowQuantiles(m, probs = 0.25, na.rm = TRUE)
  q3  <- matrixStats::rowQuantiles(m, probs = 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 3 * iqr
  upper <- q3 + 3 * iqr
  outlier_counts <- rowSums(m < lower | m > upper, na.rm = TRUE)
  
  # range and variance
  row_range <- matrixStats::rowMaxs(m, na.rm = TRUE) - matrixStats::rowMins(m, na.rm = TRUE)
  row_var   <- matrixStats::rowVars(m, na.rm = TRUE)
  
  # skewness / kurtosis
  skewness <- apply(m, 1, function(v) if (length(na.omit(v))>=3) e1071::skewness(v, type=2, na.rm = T) else NA)
  kurtosis <- apply(m, 1, function(v) if (length(na.omit(v))>=3) e1071::kurtosis(v, type=2, na.rm = T) else NA)
  
  data.frame(
    feature   = rownames(m),
    mean      = rowMeans(m, na.rm = TRUE),
    median    = matrixStats::rowMedians(m, na.rm = TRUE),
    range     = row_range,
    variance  = row_var,
    sd        = sqrt(row_var),
    min       = matrixStats::rowMins(m, na.rm = TRUE),
    max       = matrixStats::rowMaxs(m, na.rm = TRUE),
    skewness  = skewness,
    kurtosis  = kurtosis,
    n_na      = rowSums(is.na(m)),
    prop_na   = rowMeans(is.na(m)),
    q1        = q1,
    q3        = q3,
    iqr       = iqr,
    fence_lower = lower,
    fence_upper = upper,
    n_outlier   = outlier_counts,
    stringsAsFactors = FALSE
  )
}



desc_counts <- function(df_counts, rmse_ref = c("mean", "zero")) {
  # df_counts: samples x celltypes (+ optional ID column)
  # rmse_ref: whether RMSE is relative to "mean" (== sd) or "zero"
  rmse_ref <- match.arg(rmse_ref)
  
  x <- df_counts
  if ("ID" %in% names(x)) x <- x[ , setdiff(names(x), "ID"), drop = FALSE]
  x <- as.data.frame(x)
  
  tibble::tibble(
    cell_type = colnames(x),
    mean      = sapply(x, function(v) mean(v, na.rm = TRUE)),
    sd        = sapply(x, function(v) sd(v, na.rm = TRUE)),
    min       = sapply(x, function(v) min(v, na.rm = TRUE)),
    max       = sapply(x, function(v) max(v, na.rm = TRUE)),
    n_na      = sapply(x, function(v) sum(is.na(v))),
    prop_na   = sapply(x, function(v) mean(is.na(v))) 
  )
}


safe_as_matrix <- function(obj) {
  # Accepts: matrix, data.frame, meffil QC/normalized list, or list with $betas
  if (is.list(obj) && "betas" %in% names(obj)) return(as.matrix(obj$betas))
  if (is.matrix(obj)) return(obj)
  if (is.data.frame(obj)) return(as.matrix(obj))
  stop("Cannot coerce 'meth' into a beta matrix. Supply a matrix/data.frame or list with $betas.")
}

# One combined boxplot per cell-count set (samples x cell types [+ optional ID])
# Saves: boxplot_<name>_all.png
plot_cellcount_boxplot <- function(df_counts, outdir, name = "cellcounts", clamp_0_100 = TRUE) {
  if (is.null(df_counts)) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("Please install 'ggplot2' and 'tidyr'.")
  }
  x <- df_counts
  if ("ID" %in% names(x)) x$ID <- NULL
  num_cols <- vapply(x, is.numeric, logical(1))
  x <- x[, num_cols, drop = FALSE]
  if (ncol(x) == 0) {
    message("[INFO] No numeric cell-type columns in ", name, " — skipping plot.")
    return(invisible(NULL))
  }
  
  long <- data.frame(sample = seq_len(nrow(x)), x, check.names = FALSE)
  long <- tidyr::pivot_longer(long, -sample, names_to = "cell_type", values_to = "percent")
  
  p <- ggplot2::ggplot(long, ggplot2::aes(x = cell_type, y = percent)) +
    ggplot2::geom_boxplot(outlier.shape = 16, linewidth = 0.3) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Cell counts: ", name),
      x = NULL, y = "Estimated cell proportion (%)"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  
  if (clamp_0_100) {
    p <- p + ggplot2::scale_y_continuous(limits = c(0, 100))
  }
  
  # Dynamic height to avoid squishing
  h <- max(3, 0.35 * length(unique(long$cell_type)))
  ggplot2::ggsave(
    filename = file.path(outdir, paste0("boxplot_", name, "_all.png")),
    plot = p, width = 6, height = h, dpi = 300
  )
  p
}

# correlations
compute_pairwise_corr <- function(df1, df2, expra, exprb, id_col = "ID", digits = 2, verbose = TRUE) {
  stopifnot(id_col %in% names(df1), id_col %in% names(df2))

  # 1) Force suffixes so selectors are reliable
  add_suffix <- function(df, suf) {
    nms <- names(df); keep <- nms == id_col
    names(df)[!keep] <- paste0(nms[!keep], suf)
    df
  }
  dfa <- add_suffix(df1, ".a")
  dfb <- add_suffix(df2, ".b")

  # 2) Align
  m <- merge(dfa, dfb, by = id_col)

  # 3) Synonyms (system-aware)
  with_synonyms <- function(cols, system) {
    # Always de-dup at the end
    expand <- function(x) unique(x)

    if (system == "unilife") {
      # sometimes stored with _cb suffix
      base <- c("B","CD4T","CD8T","Gran","Mono","NK")
      add  <- intersect(cols, base)
      return(expand(c(cols, paste0(add, "_cb"))))
    }

    if (system %in% c("epidish_blood","hepidish_saliva")) {
      mapped <- unlist(lapply(cols, function(cn) {
        switch(cn,
          # common hepidish/epidish differences
          "CD8T"   = c("CD8T","CD8t"),
          "Eosino" = c("Eosino","Eos"),
          "Neutro" = c("Neutro","Neu"),
          "B"      = c("B","Bcell","B_cell"),
          "Bcell"  = c("Bcell","B"),
          # pass-through
          cn
        )
      }), use.names = FALSE)
      return(expand(mapped))
    }

    cols
  }

  pick_present <- function(side_cols, suffix, system) {
    wanted   <- with_synonyms(side_cols, system)
    present  <- paste0(wanted, suffix)
    present[present %in% colnames(m)]
  }

  safe_sum <- function(df, cols) {
    if (length(cols) == 0) return(rep(NA_real_, nrow(df)))
    rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
  }

  build_rows <- function(mapping, comp_a, comp_b) {
    rows <- lapply(names(mapping), function(ct) {
      want_a <- mapping[[ct]]$a
      want_b <- mapping[[ct]]$b
      have_a <- pick_present(want_a, ".a", comp_a)
      have_b <- pick_present(want_b, ".b", comp_b)

      if (length(have_a) == 0 || length(have_b) == 0) {
        if (verbose) {
          missing_side <- if (length(have_a) == 0 && length(have_b) == 0) "both"
          else if (length(have_a) == 0) "A" else "B"
          message(sprintf(
            "[compute_pairwise_corr] Skipping '%s' for %s vs %s: missing required columns on side %s.\n  needed A=%s | B=%s\n  have   A=%s | B=%s",
            ct, comp_a, comp_b, missing_side,
            paste(want_a, collapse = ", "),
            paste(want_b, collapse = ", "),
            paste(sub("\\.a$", "", have_a), collapse = ", "),
            paste(sub("\\.b$", "", have_b), collapse = ", ")
          ))
        }
        return(NULL)
      }

      a_vals <- safe_sum(m, have_a)
      b_vals <- safe_sum(m, have_b)
      ok <- is.finite(a_vals) & is.finite(b_vals)
      n_ok <- sum(ok)
      if (n_ok == 0) return(NULL)

      r_val <- suppressWarnings(cor(a_vals[ok], b_vals[ok]))
      rmse  <- sqrt(mean((a_vals[ok] - b_vals[ok])^2))

      data.frame(
        cell_type = ct,
        r    = round(r_val, digits),
        rmse = round(rmse,  digits),
        n    = n_ok,
        comp_a = comp_a,
        comp_b = comp_b,
        expr_a = paste(want_a, collapse = " + "),
        expr_b = paste(want_b, collapse = " + "),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }

  # 4) Mappings ----------------------------------------------------------
  if (expra == "meffil" & exprb == "unilife") {
    mapping <- list(
      Bcell = list(a = c("Bcell"),                    b = c("aBmem","aBnv","B")),
      CD4T  = list(a = c("CD4T"),                     b = c("CD4T","aCD4Tnv","aCD4Tmem","aTreg")),
      CD8T  = list(a = c("CD8T"),                     b = c("CD8T","aCD8Tnv","aCD8Tmem")),
      Gran  = list(a = c("Gran"),                     b = c("aEos","aNeu","aBaso","Gran")),
      Mono  = list(a = c("Mono"),                     b = c("aMono","Mono")),
      NK    = list(a = c("NK"),                       b = c("aNK","NK")),
      nRBC  = list(a = c("nRBC"),                       b = c("nRBC"))
    )
    return(build_rows(mapping, "meffil", "unilife"))
  }

  if (expra == "unilife" & exprb == "epidish_blood") {
    mapping <- list(
      Mono     = list(a = c("aMono","Mono"),          b = c("Mono")),
      Gran     = list(a = c("aEos","aNeu","aBaso","Gran"), b = c("Eos","Neu","Baso")),
      NK       = list(a = c("aNK","NK"),              b = c("NK")),
      Bnv      = list(a = c("B","aBnv"),              b = c("Bnv")),
      Bmem     = list(a = c("aBmem"),                 b = c("Bmem")),
      CD4Tnv   = list(a = c("aCD4Tnv","CD4T"),        b = c("CD4Tnv")),
      Treg     = list(a = c("aTreg"),                 b = c("Treg")),
      CD4Tmem  = list(a = c("aCD4Tmem"),              b = c("CD4Tmem")),
      CD8Tnv   = list(a = c("CD8T","aCD8Tnv"),        b = c("CD8Tnv")),
      CD8Tmem  = list(a = c("aCD8Tmem"),              b = c("CD8Tmem"))
    )
    return(build_rows(mapping, "unilife", "epidish_blood"))
  }

  if (expra == "meffil" & exprb == "epidish_blood") {
    mapping <- list(
      Bcell = list(a = c("Bcell"),                    b = c("Bnv","Bmem")),
      CD4T  = list(a = c("CD4T"),                     b = c("CD4Tnv","CD4Tmem","Treg")),
      CD8T  = list(a = c("CD8T"),                     b = c("CD8Tnv","CD8Tmem")),
      Gran  = list(a = c("Gran"),                     b = c("Eos","Neu","Baso")),
      Mono  = list(a = c("Mono"),                     b = c("Mono")),
      NK    = list(a = c("NK"),                       b = c("NK"))
    )
    return(build_rows(mapping, "meffil", "epidish_blood"))
  }

  # -------- NEW: uniLIFE vs hepiDISH (saliva/buccal) --------
  if (expra == "unilife" & exprb == "hepidish_saliva") {
    mapping <- list(
      Mono   = list(a = c("aMono","Mono"),                      b = c("Mono")),
      Gran   = list(a = c("aEos","aNeu","aBaso","Gran"),        b = c("Eosino","Neutro")), # Eos + Neu
      NK     = list(a = c("aNK","NK"),                          b = c("NK")),
      B      = list(a = c("B","aBnv","aBmem"),                  b = c("B")),              # B + aBnv + Bmem = B
      CD4T   = list(a = c("CD4T","aCD4Tnv","aCD4Tmem","aTreg"), b = c("CD4T")),           # CD4T composite = CD4T
      CD8T   = list(a = c("CD8T","aCD8Tnv","aCD8Tmem"),         b = c("CD8T"))            # CD8T composite = CD8T
    )
    return(build_rows(mapping, "unilife", "hepidish_saliva"))
  }

  # -------- NEW: meffil vs hepiDISH (saliva/buccal) --------
  if (expra == "meffil" & exprb == "hepidish_saliva") {
    mapping <- list(
      CD4T   = list(a = c("CD4T"),            b = c("CD4T")),
      CD8T   = list(a = c("CD8T"),            b = c("CD8T")),
      Gran   = list(a = c("Gran"),            b = c("Eosino","Neutro")),  # Eos + Neu
      Mono   = list(a = c("Mono"),            b = c("Mono")),
      NK     = list(a = c("NK"),              b = c("NK")),
      Buccal = list(a = c("Buccal"),          b = c("Epi","Fib")),        # report Epi+Fib as 'Buccal'
      Bcell  = list(a = c("Bcell"),           b = c("B"))                 # Bcell = B
    )
    return(build_rows(mapping, "meffil", "hepidish_saliva"))
  }

  # Optional: unilife vs unilife (self-consistency)
  if (expra == "unilife" & exprb == "unilife") {
    mapping <- list(
      Bcell = list(a = c("aBmem","aBnv","B"),         b = c("aBmem","aBnv","B")),
      CD4T  = list(a = c("CD4T","aCD4Tnv","aCD4Tmem","aTreg"),
                   b = c("CD4T","aCD4Tnv","aCD4Tmem","aTreg")),
      CD8T  = list(a = c("CD8T","aCD8Tnv","aCD8Tmem"),
                   b = c("CD8T","aCD8Tnv","aCD8Tmem")),
      Gran  = list(a = c("aEos","aNeu","aBaso","Gran"),
                   b = c("aEos","aNeu","aBaso","Gran")),
      Mono  = list(a = c("aMono","Mono"),             b = c("aMono","Mono")),
      NK    = list(a = c("aNK","NK"),                 b = c("aNK","NK"))
    )
    return(build_rows(mapping, "unilife", "unilife"))
  }

  stop(sprintf("Unsupported expra/exprb combination: %s vs %s", expra, exprb))
}
