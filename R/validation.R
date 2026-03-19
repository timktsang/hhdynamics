# Input validation for household_dynamics()

validate_inputs <- function(input, inf_factor, sus_factor, SI, n_iteration, burnin, thinning, with_rm) {
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

  # Check index cases are infected
  index_rows <- input$member == 0
  if (any(index_rows)) {
    if (any(input$inf[index_rows] != 1, na.rm = TRUE)) {
      stop("All index cases (member == 0) must be infected (inf == 1).", call. = FALSE)
    }
  }

  # Check formula arguments
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

  # Check SI
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
