
#' template.match
#'
#' Compute correlation between one feature vector and a template
#'
#' @param x         numeric vector — one feature, all samples
#' @param template  numeric vector — same length as x, the reference pattern
#' @param dist.name character(1) — "pearson", "spearman", or "kendall"
#' @return numeric(3): c(estimate, statistic, p.value)
#' @keywords internal
template.match <- function(x, template, dist.name) {
  k <- cor.test(x, template, method = dist.name, exact = FALSE)
  c(k$estimate, k$statistic, k$p.value)
}


#' run_pattern_search_from_gct
#' norm_mat rows and metadata_df rows must already be aligned (same samples,
#' same order) — PollyCommonR::sample_intensity_matrix() in the server
#' guarantees this before this function is called.
#'
#' @param norm_mat     matrix, samples x features (already transposed from @mat)
#' @param metadata_df  data.frame, samples x metadata cols (filtered @cdesc)
#' @param cohort_col   character(1) — which column holds group labels
#' @param pattern      character(1) — e.g. "1-2-3", NULL auto-generates
#' @param dist.name    character(1) — "pearson", "spearman", or "kendall"
#' @return list(cor_mat, pattern, method, groups) or stops with message
#' @keywords internal
run_pattern_search_from_gct <- function(norm_mat,
                                         metadata_df,
                                         cohort_col,
                                         pattern   = NULL,
                                         dist.name = "pearson") {

  if (is.null(norm_mat) || nrow(norm_mat) == 0)
    stop("norm_mat is NULL or empty.")
  if (is.null(metadata_df) || nrow(metadata_df) == 0)
    stop("metadata_df is NULL or empty.")
  if (!cohort_col %in% colnames(metadata_df))
    stop(paste0(
      "Cohort column '", cohort_col, "' not found in metadata. ",
      "Available columns: ", paste(colnames(metadata_df), collapse = ", ")
    ))

  # Build cls factor.
  # Level order = unique() order of the filtered metadata — which reflects
  # the order the user selected cohorts in ptsr01_cohort_selected.
  # This is what maps pattern position 1 → first cohort, 2 → second, etc.
  cls      <- factor(
    metadata_df[[cohort_col]],
    levels = sort(unique(metadata_df[[cohort_col]]))
  )
  names(cls) <- rownames(metadata_df)
  all.lvls   <- levels(cls)

  # Default pattern
  if (is.null(pattern) || nchar(trimws(pattern)) == 0)
    pattern <- paste(seq_along(all.lvls), collapse = "-")

  # Parse and validate
  templ <- suppressWarnings(
    as.numeric(trimws(strsplit(pattern, "-", fixed = TRUE)[[1]]))
  )

  if (any(is.na(templ)))
    stop(paste0(
      "Pattern '", pattern, "' contains non-numeric values. ",
      "Use integers separated by '-', e.g. '1-2-3'."
    ))
  if (all(templ == templ[1]))
    stop(
      "Cannot correlate on a constant pattern (all values identical). ",
      "Each group must have a different pattern value."
    )
  if (length(templ) != length(all.lvls))
    stop(paste0(
      "Pattern length mismatch: pattern '", pattern, "' has ", length(templ),
      " value(s) but '", cohort_col, "' has ", length(all.lvls),
      " group(s): ", paste(all.lvls, collapse = ", "),
      ".\nPattern must have exactly one value per group."
    ))

  # Expand group-level template to sample-level vector
  new.template <- numeric(length(cls))
  for (i in seq_along(all.lvls))
    new.template[cls == all.lvls[i]] <- templ[i]

  # Align norm_mat rows to metadata_df (@cdesc) row order.
  common_samples <- names(cls)[names(cls) %in% rownames(norm_mat)]
  if (length(common_samples) == 0)
    stop(
      "No common samples between norm_mat and metadata_df. ",
      "Check that rownames(norm_mat) match rownames(metadata_df)."
    )
  norm_mat     <- norm_mat[common_samples, , drop = FALSE]
  new.template <- new.template[match(common_samples, names(cls))]

  # Correlate
  is.partial <- grepl("^p-|^partial_", dist.name)

  if (is.partial) {
    if (!requireNamespace("ppcor", quietly = TRUE))
      stop("Package 'ppcor' is required for partial correlation. ",
           "Run: install.packages('ppcor')")
    method     <- sub("^(p-|partial_)", "", dist.name)
    method     <- match.arg(method, c("pearson", "spearman"))
    feat.names <- colnames(norm_mat)
    cor.res    <- matrix(NA_real_, nrow = ncol(norm_mat), ncol = 3,
                         dimnames = list(feat.names,
                                         c("correlation", "t-stat", "p-value")))
    for (j in seq_along(feat.names)) {
      x_j      <- norm_mat[, j]
      z_j      <- norm_mat[, -j, drop = FALSE]
      complete <- complete.cases(x_j, new.template, z_j)
      if (sum(complete) < 5) next
      ans <- tryCatch(
        ppcor::pcor.test(x_j[complete], new.template[complete],
                         z_j[complete, ], method = method),
        error = function(e) NULL
      )
      if (!is.null(ans))
        cor.res[j, ] <- c(ans$estimate, ans$statistic, ans$p.value)
    }
  } else {
    tmp     <- apply(norm_mat, 2, template.match, new.template, dist.name)
    cor.res <- t(tmp)
  }

  # FDR, sort, round
  rownames(cor.res) <- colnames(norm_mat)
  cor.res <- cor.res[!is.na(cor.res[, 3]), , drop = FALSE]
  cor.res <- cbind(cor.res, FDR = p.adjust(cor.res[, 3], method = "fdr"))
  colnames(cor.res) <- c("correlation", "t-stat", "p-value", "FDR")
  sig.mat <- signif(cor.res[order(cor.res[, 3]), ], 5)

  return(list(
    cor_mat = sig.mat,
    pattern = pattern,
    method  = dist.name,
    groups  = all.lvls
  ))
}


#' run_feature_correlation_from_gct
#' @param norm_mat      matrix, samples x features (filter_pattern_norm_mat)
#'                      rows = samples, cols = features
#'                      Must already be filtered and aligned by the server
#'                      (same matrix used by run_pattern_search_from_gct)
#' @param query_feature character(1) — name of the reference metabolite
#'                      Must be a colname of norm_mat
#' @param dist.name     character(1) — "pearson" | "spearman" | "kendall"
#'                      Prefix "p-" or "partial_" for partial correlation
#' @return list(cor_mat, pattern, method, groups) — same structure as
#'         run_pattern_search_from_gct so the run button handler is identical
#' @keywords internal
run_feature_correlation_from_gct <- function(norm_mat,
                                              query_feature,
                                              dist.name = "pearson") {

  # ── Validate ──────────────────────────────────────────────────────────────
  if (is.null(norm_mat) || nrow(norm_mat) == 0)
    stop("norm_mat is NULL or empty.")
  if (!query_feature %in% colnames(norm_mat))
    stop(paste0(
      "Metabolite '", query_feature, "' not found in processed data. ",
      "It may have been filtered out during preprocessing. ",
      "Available features: ", paste(head(colnames(norm_mat), 5), collapse = ", "), "..."
    ))

  # ── Extract template vector ───────────────────────────────────────────────
  target.vec <- norm_mat[, query_feature]

  if (all(is.na(target.vec)))
    stop(paste0("Metabolite '", query_feature, "' has all NA values in the selected samples."))
  if (sd(target.vec, na.rm = TRUE) == 0)
    stop(paste0("Metabolite '", query_feature, "' has zero variance — cannot correlate."))

  # ── Correlate all features against the query feature ─────────────────────
  is.partial <- grepl("^p-|^partial_", dist.name)

  if (is.partial) {
    if (!requireNamespace("ppcor", quietly = TRUE))
      stop("Package 'ppcor' is required for partial correlation. ",
           "Run: install.packages('ppcor')")
    method   <- sub("^(p-|partial_)", "", dist.name)
    method   <- match.arg(method, c("pearson", "spearman"))

    # Exclude query feature from the search space for partial correlation
    features <- setdiff(colnames(norm_mat), query_feature)
    cor.res  <- matrix(NA_real_, nrow = length(features), ncol = 3,
                       dimnames = list(features,
                                       c("correlation", "t-stat", "p-value")))

    for (i in seq_along(features)) {
      test_vec  <- norm_mat[, features[i]]
      cond_vars <- norm_mat[, setdiff(features, features[i]), drop = FALSE]
      complete  <- complete.cases(test_vec, target.vec, cond_vars)
      if (sum(complete) < 5) next
      ans <- tryCatch(
        ppcor::pcor.test(
          x = test_vec[complete],
          y = target.vec[complete],
          z = cond_vars[complete, , drop = FALSE],
          method = method
        ),
        error = function(e) NULL
      )
      if (!is.null(ans))
        cor.res[i, ] <- c(ans$estimate, ans$statistic, ans$p.value)
    }

  } else {
    # Ordinary correlation — vectorised, fast.
    tmp     <- apply(norm_mat, 2, template.match, target.vec, dist.name)
    cor.res <- t(tmp)
  }

  # ── FDR, sort, round — identical to run_pattern_search_from_gct ──────────
  cor.res <- cor.res[!is.na(cor.res[, 3]), , drop = FALSE]
  cor.res <- cbind(cor.res, FDR = p.adjust(cor.res[, 3], method = "fdr"))
  colnames(cor.res) <- c("correlation", "t-stat", "p-value", "FDR")
  sig.mat <- signif(cor.res[order(cor.res[, 3]), ], 5)

  # ── Return same list structure as run_pattern_search_from_gct ─────────────
  # pattern = query feature name (for download filename, report, etc.)
  # groups  = NULL (no group structure used in this mode)
  return(list(
    cor_mat = sig.mat,
    pattern = query_feature,
    method  = dist.name,
    groups  = NULL
  ))
}


#' compute_pattern_search
#'
#' Shiny-facing wrapper called from observeEvent in the server above.
#' Reads @mat and @cdesc from norm_gct_object 
#' In normal server flow this is NOT called directly — the observeEvent above
#' handles the full flow including filter_metadata_by_cohorts and
#' sample_intensity_matrix before calling run_pattern_search_from_gct.
#' This function is provided for cases where you want to call the algorithm
#' directly without the full filter flow (e.g. scripted batch runs).
#'
#' @param norm_gct_object  the GCT object (metab_reactive$norm_gct_object)
#' @param cohort_col       character(1) — which @cdesc column holds groups
#' @param pattern          character(1) — e.g. "1-2-3"
#' @param dist.name        character(1) — correlation method
#' @return list(cor_mat, pattern, method, groups) or NULL on error
#' @keywords internal
compute_pattern_search <- function(norm_gct_object,
                                    cohort_col,
                                    pattern   = NULL,
                                    dist.name = "pearson") {

  if (is.null(norm_gct_object)) {
    message("compute_pattern_search: norm_gct_object is NULL")
    return(NULL)
  }

  norm_mat_df     <- as.data.frame(
    norm_gct_object@mat,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  norm_mat_df$id  <- rownames(norm_mat_df)
  metadata_df     <- as.data.frame(
    norm_gct_object@cdesc,
    stringsAsFactors = FALSE, check.names = FALSE
  )

  norm_mat_aligned <- tryCatch(
    PollyCommonR::sample_intensity_matrix(
      raw_intensity_df = norm_mat_df,
      metadata_df      = metadata_df,
      rownames_col     = "id"
    ),
    error = function(e) {
      message("sample_intensity_matrix failed, using manual transpose: ", e$message)
      t(norm_gct_object@mat)
    }
  )

  tryCatch(
    run_pattern_search_from_gct(
      norm_mat    = norm_mat_aligned,
      metadata_df = metadata_df,
      cohort_col  = cohort_col,
      pattern     = pattern,
      dist.name   = dist.name
    ),
    error = function(e) {
      message("Pattern search error: ", e$message)
      return(NULL)
    }
  )
}


#' create_pattern_search_markdown
#'
#' Builds the markdown string for the parameters report entry.
#' Called in the run button handler immediately after results are stored.
#'
#' @param mode        "manual" or "feature"
#' @param cohort_col  character — selected cohort column name
#' @param cohorts     character vector — selected cohort values
#' @param pattern     character — "1-2-3" for manual, metabolite name for feature
#' @param dist_method character — correlation method used
#' @param n_features  integer — number of features in the result
#' @return character(1) — markdown string
#' @keywords internal
create_pattern_search_markdown <- function(mode        = NULL,
                                            cohort_col  = NULL,
                                            cohorts     = NULL,
                                            pattern     = NULL,
                                            dist_method = "pearson",
                                            n_features  = NULL) {

  mode_label <- if (!is.null(mode) && mode == "feature") {
    "Metabolite of Interest"
  } else {
    "Manual Pattern"
  }

  pattern_label <- if (!is.null(mode) && mode == "feature") {
    paste0("Query metabolite: ", pattern)
  } else {
    paste0("Pattern: ", pattern)
  }

  paste0(
    "# Pattern Search",
    "\n * Mode: ",               mode_label,
    "\n * ", pattern_label,
    "\n * Cohort column: ",      ifelse(!is.null(cohort_col), cohort_col, "N/A"),
    "\n * Cohorts: ",            ifelse(!is.null(cohorts),
                                        paste(cohorts, collapse = ", "), "N/A"),
    "\n * Correlation method: ", dist_method,
    "\n * Features in result: ", ifelse(!is.null(n_features), n_features, "N/A"),
    "\n"
  )
}


#' build_correlation_bar_chart
#'
#' Create a horizontal bar chart for correlation results
#'
#' @param cor_mat correlation results matrix with columns "correlation", "p-value", "FDR"
#' @param title chart title string
#' @return plotly object
#' @details
#' Takes the top 25 results by p-value and re-sorts by correlation magnitude.
#' Positive correlations shown in red, negative in cyan.
#' @keywords internal
build_correlation_bar_chart <- function(cor_mat, title = "Pattern Search Results") {

  if (is.null(cor_mat) || nrow(cor_mat) == 0)
    stop("cor_mat is NULL or empty.")

  # Top 25 by p-value (cor_mat already sorted by p-value from the algorithm)
  plot_mat <- if (nrow(cor_mat) > 25) cor_mat[1:25, , drop = FALSE] else cor_mat

  # Re-sort by correlation value — mirrors PlotCorr() logic
  ord <- order(plot_mat[, "correlation"])
  if (all(plot_mat[, "correlation"] < 0)) ord <- rev(ord)
  plot_mat <- plot_mat[ord, , drop = FALSE]

  features   <- rownames(plot_mat)
  corr_vals  <- as.numeric(plot_mat[, "correlation"])
  pval_vals  <- as.numeric(plot_mat[, "p-value"])
  bar_colors <- ifelse(corr_vals >= 0, "#E64B35", "#4DBBD5")

  hover_text <- paste0(
    "<b>", features, "</b><br>",
    "Correlation: ", round(corr_vals, 4), "<br>",
    "p-value: ",     signif(pval_vals, 3), "<br>"
  )

  plotly::plot_ly(
    x = corr_vals,
    y = factor(features, levels = features),
    type = "bar", orientation = "h",
    marker = list(color = bar_colors),
    hoverinfo = "text", text = hover_text
  ) %>%
    plotly::layout(
      title = list(text = title, font = list(size = 14)),
      xaxis = list(title = "Correlation coefficient", range = c(-1, 1),
                   zeroline = TRUE, zerolinecolor = "#888888", zerolinewidth = 1.5),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 10, r = 30, t = 50, b = 50),
      bargap = 0.3,
      plot_bgcolor = "white", paper_bgcolor = "white"
    )
}

