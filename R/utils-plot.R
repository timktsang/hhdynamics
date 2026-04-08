# Internal helper functions for plotting and tables

# Effective sample size using initial positive sequence estimator (Geyer 1992)
.effective_sample_size <- function(x) {
  n <- length(x)
  if (n < 4) return(n)
  if (stats::var(x) == 0) return(1)
  acf_vals <- stats::acf(x, lag.max = n - 1, plot = FALSE)$acf[, , 1]
  # Sum consecutive pairs of autocorrelations starting at lag 1
  # (acf_vals[1] is lag 0 = 1.0, so skip it)
  tau <- 1
  k <- 2  # start at lag 1
  while (k + 1 <= length(acf_vals)) {
    pair_sum <- acf_vals[k] + acf_vals[k + 1]
    if (is.na(pair_sum) || pair_sum < 0) break
    tau <- tau + 2 * pair_sum
    k <- k + 2
  }
  ess <- n / tau
  max(1, min(ess, n))
}

# Transform parameter samples for display
.transform_param <- function(x, transform) {
  switch(transform,
    "prob" = 1 - exp(-x),
    "exp"  = exp(x),
    x  # "none" or anything else
  )
}

# Wilson score confidence interval for a binomial proportion
.wilson_ci <- function(k, n, alpha = 0.05) {
  if (n == 0) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(1 - alpha / 2)
  p_hat <- k / n
  denom <- 1 + z^2 / n
  centre <- (p_hat + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * n)) / n) / denom
  c(lower = max(0, centre - margin), upper = min(1, centre + margin))
}

# Build structured row data for forest plot
.build_forest_rows <- function(fit, probs, labels = NULL) {
  rows <- list()

  for (type in c("sus", "inf")) {
    formula <- if (type == "inf") fit$formula_inf else fit$formula_sus
    n_type <- if (type == "inf") fit$n_inf else fit$n_sus
    if (is.null(formula) || n_type == 0) next

    section_label <- if (type == "sus") "Factors Affecting Susceptibility:"
                     else "Factors Affecting Infectiousness:"
    rows[[length(rows) + 1]] <- list(type = "section", label = section_label)

    # Param offset in samples matrix: inf at 5:(4+n_inf), sus at (5+n_inf):...
    param_offset <- if (type == "inf") 4 else 4 + fit$n_inf

    # Reconstruct design matrix structure for term mapping
    mf <- stats::model.frame(formula, fit$input_data, na.action = stats::na.pass)
    dm <- model.matrix(attr(mf, "terms"), mf)
    assign_vec <- attr(dm, "assign")[-1]  # drop intercept
    term_labels <- attr(stats::terms(formula), "term.labels")

    col_counter <- 0
    for (term_idx in seq_along(term_labels)) {
      var_name <- term_labels[term_idx]
      n_cols <- sum(assign_vec == term_idx)
      var_data <- fit$input_data[[var_name]]

      # Custom labels
      var_display <- var_name
      custom_levels <- NULL
      if (!is.null(labels) && var_name %in% names(labels)) {
        lab <- labels[[var_name]]
        if (!is.null(lab$name)) var_display <- lab$name
        if (!is.null(lab$levels)) custom_levels <- lab$levels
      }

      is_fac <- is.factor(var_data) || is.character(var_data)

      if (is_fac) {
        levs <- levels(factor(var_data))
        level_labels <- if (!is.null(custom_levels)) custom_levels else levs

        # Variable header row
        rows[[length(rows) + 1]] <- list(type = "header", label = var_display)
        # Reference level row
        rows[[length(rows) + 1]] <- list(type = "ref", label = level_labels[1])

        # Non-reference level rows with estimates
        for (j in seq_len(n_cols)) {
          col_counter <- col_counter + 1
          pi <- param_offset + col_counter
          x <- fit$samples[, pi]
          est <- exp(mean(x))
          lo <- exp(stats::quantile(x, probs[1], names = FALSE))
          hi <- exp(stats::quantile(x, probs[2], names = FALSE))
          lbl <- if (j + 1 <= length(level_labels)) level_labels[j + 1] else levs[j + 1]

          rows[[length(rows) + 1]] <- list(
            type = "estimate", label = lbl,
            estimate = est, lower = lo, upper = hi, is_level = TRUE
          )
        }
      } else {
        # Continuous covariate: single estimate row, no header/ref
        col_counter <- col_counter + 1
        pi <- param_offset + col_counter
        x <- fit$samples[, pi]
        est <- exp(mean(x))
        lo <- exp(stats::quantile(x, probs[1], names = FALSE))
        hi <- exp(stats::quantile(x, probs[2], names = FALSE))

        rows[[length(rows) + 1]] <- list(
          type = "estimate", label = var_display,
          estimate = est, lower = lo, upper = hi, is_level = FALSE
        )
      }
    }
  }

  rows
}

# Get parameter indices to plot (skip fixed params)
.get_plot_params <- function(fit, params = NULL) {
  all_names <- fit$param_names
  # Skip fixed params by default
  skip <- "size_param"
  if (fit$with_rm == 0) skip <- c(skip, "re_sd")
  default_names <- setdiff(all_names, skip)

  if (!is.null(params)) {
    bad <- setdiff(params, all_names)
    if (length(bad) > 0) {
      stop("Unknown parameter names: ", paste(bad, collapse = ", "),
           ". Available: ", paste(all_names, collapse = ", "), call. = FALSE)
    }
    return(params)
  }
  default_names
}
