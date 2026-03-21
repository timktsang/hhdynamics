# Input validation for household_dynamics()

validate_inputs <- function(input, inf_factor, sus_factor, SI, n_iteration, burnin, thinning, with_rm, estimate_SI = FALSE) {
  # Check input is a data frame
  if (!is.data.frame(input)) {
    stop("'input' must be a data frame.", call. = FALSE)
  }

  # Check required columns
  required_cols <- c("hhID", "member", "size", "end", "inf", "onset")
  missing_cols <- setdiff(required_cols, names(input))
  if (length(missing_cols) > 0) {
    stop(sprintf("'input' is missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Check no missing values in required columns (onset checked separately below)
  for (col in setdiff(required_cols, "onset")) {
    if (any(is.na(input[[col]]))) {
      n_na <- sum(is.na(input[[col]]))
      stop(sprintf("Column '%s' has %d missing value(s). This version requires complete data (no NAs in required columns).",
                    col, n_na), call. = FALSE)
    }
  }

  # Check index cases are infected
  index_rows <- input$member == 0
  if (any(index_rows)) {
    if (any(input$inf[index_rows] != 1)) {
      stop("All index cases (member == 0) must be infected (inf == 1).", call. = FALSE)
    }
  }

  # Check infected contacts have non-missing onset times
  infected_contacts <- input$member != 0 & input$inf == 1
  if (any(infected_contacts)) {
    bad_onset <- infected_contacts & (is.na(input$onset) | input$onset < 0)
    if (any(bad_onset)) {
      n_bad <- sum(bad_onset)
      stop(sprintf("%d infected contact(s) have missing or negative onset times. This version requires known onset for all infected individuals.",
                    n_bad), call. = FALSE)
    }
  }

  # Check formula arguments (must come before any all.vars/terms calls)
  if (!is.null(inf_factor)) {
    if (!inherits(inf_factor, "formula")) {
      stop("'inf_factor' must be a formula (e.g. ~sex) or NULL.", call. = FALSE)
    }
    formula_vars <- all.vars(inf_factor)
    missing_vars <- setdiff(formula_vars, names(input))
    if (length(missing_vars) > 0) {
      stop(sprintf("Variables in inf_factor not found in input: %s", paste(missing_vars, collapse = ", ")), call. = FALSE)
    }
  }

  if (!is.null(sus_factor)) {
    if (!inherits(sus_factor, "formula")) {
      stop("'sus_factor' must be a formula (e.g. ~age) or NULL.", call. = FALSE)
    }
    formula_vars <- all.vars(sus_factor)
    missing_vars <- setdiff(formula_vars, names(input))
    if (length(missing_vars) > 0) {
      stop(sprintf("Variables in sus_factor not found in input: %s", paste(missing_vars, collapse = ", ")), call. = FALSE)
    }
  }

  # Check formula covariates — allow NAs for factors (MCMC will impute), error for continuous
  covariate_vars <- unique(c(
    if (!is.null(inf_factor)) all.vars(inf_factor),
    if (!is.null(sus_factor)) all.vars(sus_factor)
  ))
  for (v in covariate_vars) {
    if (v %in% names(input) && any(is.na(input[[v]]))) {
      if (!is.factor(input[[v]])) {
        n_na <- sum(is.na(input[[v]]))
        stop(sprintf("Continuous covariate '%s' has %d missing value(s). Only factor covariates can be imputed. Convert to factor or remove NAs.",
                      v, n_na), call. = FALSE)
      }
      n_na <- sum(is.na(input[[v]]))
      pct <- round(100 * n_na / nrow(input), 1)
      message(sprintf("Note: Covariate '%s' has %d missing value(s) (%.1f%%). These will be imputed during MCMC.",
                      v, n_na, pct))
    }
  }

  # Check for interaction terms with missing covariate data
  for (f in list(inf_factor, sus_factor)) {
    if (!is.null(f)) {
      tl <- attr(stats::terms(f), "term.labels")
      if (any(grepl(":", tl))) {
        for (v in all.vars(f)) {
          if (v %in% names(input) && any(is.na(input[[v]]))) {
            stop("Interaction terms with missing covariate data are not supported. Remove interactions or pre-impute missing values.",
                 call. = FALSE)
          }
        }
      }
    }
  }

  # Check for shared variables with missing values across formulas
  if (!is.null(inf_factor) && !is.null(sus_factor)) {
    shared_vars <- intersect(all.vars(inf_factor), all.vars(sus_factor))
    for (v in shared_vars) {
      if (v %in% names(input) && any(is.na(input[[v]]))) {
        stop(sprintf(paste0("Covariate '%s' appears in both inf_factor and sus_factor and has ",
                            "missing values. This is not supported. Use the variable in only one ",
                            "formula, or pre-impute the missing values."), v), call. = FALSE)
      }
    }
  }

  # Check SI (skip when estimating SI from data)
  if (!estimate_SI) {
    if (!is.numeric(SI)) {
      stop("'SI' must be a numeric vector.", call. = FALSE)
    }
    if (length(SI) != 14) {
      stop(sprintf("'SI' must have length 14, got %d.", length(SI)), call. = FALSE)
    }
    if (any(SI < 0)) {
      stop("'SI' must have all non-negative values.", call. = FALSE)
    }
    si_sum <- sum(SI)
    if (abs(si_sum - 1) > 0.01) {
      stop(sprintf("'SI' must sum to approximately 1, got %.4f.", si_sum), call. = FALSE)
    }
  }

  # Check MCMC parameters
  if (!is.numeric(n_iteration) || n_iteration < 1) {
    stop("'n_iteration' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(burnin) || burnin < 0) {
    stop("'burnin' must be a non-negative integer.", call. = FALSE)
  }
  if (n_iteration <= burnin) {
    stop("'n_iteration' must be greater than 'burnin'.", call. = FALSE)
  }
  if (!is.numeric(thinning) || thinning < 1) {
    stop("'thinning' must be a positive integer.", call. = FALSE)
  }

  # Report max household size for awareness
  max_hh <- max(input$size, na.rm = TRUE)
  if (max_hh > 20) {
    message(sprintf("Note: maximum household size is %d.", max_hh))
  }

  invisible(TRUE)
}
