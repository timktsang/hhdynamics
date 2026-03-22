# Table functions for hhdynamics

#' Table of transmission parameter estimates
#'
#' Returns a clean data frame of posterior parameter estimates with credible
#' intervals and effective sample sizes. Community and household parameters
#' are reported on the probability scale.
#'
#' @param fit An object of class \code{hhdynamics_fit}.
#' @param probs Numeric vector of length 2 for credible interval bounds.
#'   Default: \code{c(0.025, 0.975)} (95\% CrI).
#' @return A data frame with columns: Parameter, Mean, Median, Lower, Upper,
#'   ESS, Acceptance.
#' @examples
#' \donttest{
#' data(inputdata)
#' fit <- household_dynamics(inputdata, ~sex, ~age,
#'   n_iteration = 15000, burnin = 5000, thinning = 1)
#' table_parameters(fit)
#' }
#' @export
table_parameters <- function(fit, probs = c(0.025, 0.975)) {
  plot_names <- .get_plot_params(fit, NULL)
  idx <- match(plot_names, fit$param_names)
  n <- length(plot_names)

  out <- data.frame(
    Parameter = character(n),
    Mean = numeric(n),
    Median = numeric(n),
    Lower = numeric(n),
    Upper = numeric(n),
    ESS = integer(n),
    Acceptance = numeric(n),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n)) {
    x_raw <- fit$samples[, idx[i]]
    transform <- fit$param_transform[idx[i]]
    x_display <- .transform_param(x_raw, transform)

    display_name <- switch(plot_names[i],
      "community" = "Daily P(community infection)",
      "household" = "P(household transmission)",
      plot_names[i]
    )

    out$Parameter[i]  <- display_name
    out$Mean[i]       <- mean(x_display)
    out$Median[i]     <- stats::median(x_display)
    out$Lower[i]      <- stats::quantile(x_display, probs[1])
    out$Upper[i]      <- stats::quantile(x_display, probs[2])
    out$ESS[i]        <- round(.effective_sample_size(x_raw))
    out$Acceptance[i]  <- fit$acceptance[idx[i]]
  }

  rownames(out) <- NULL
  out
}

#' Table of covariate effects (odds ratios)
#'
#' Returns a data frame of covariate effects on infectivity and susceptibility,
#' with exponentiated estimates interpretable as relative risks.
#'
#' @param fit An object of class \code{hhdynamics_fit}.
#' @param probs Numeric vector of length 2 for credible interval bounds.
#' @return A data frame with columns: Covariate, Type, Estimate, Lower, Upper,
#'   exp_Estimate, exp_Lower, exp_Upper. Returns an empty data frame if no
#'   covariates were used.
#' @examples
#' \donttest{
#' data(inputdata)
#' fit <- household_dynamics(inputdata, ~sex, ~age,
#'   n_iteration = 15000, burnin = 5000, thinning = 1)
#' table_covariates(fit)
#' }
#' @export
table_covariates <- function(fit, probs = c(0.025, 0.975)) {
  n_cov <- fit$n_inf + fit$n_sus
  if (n_cov == 0) {
    return(data.frame(
      Covariate = character(0), Type = character(0),
      Estimate = numeric(0), Lower = numeric(0), Upper = numeric(0),
      exp_Estimate = numeric(0), exp_Lower = numeric(0), exp_Upper = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  cov_idx <- 5:(4 + n_cov)
  cov_names <- fit$param_names[cov_idx]
  n <- length(cov_names)

  types <- character(n)
  if (fit$n_inf > 0) types[seq_len(fit$n_inf)] <- "Infectivity"
  if (fit$n_sus > 0) types[(fit$n_inf + 1):n] <- "Susceptibility"

  out <- data.frame(
    Covariate = cov_names,
    Type = types,
    Estimate = numeric(n),
    Lower = numeric(n),
    Upper = numeric(n),
    exp_Estimate = numeric(n),
    exp_Lower = numeric(n),
    exp_Upper = numeric(n),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n)) {
    x <- fit$samples[, cov_idx[i]]
    out$Estimate[i] <- mean(x)
    out$Lower[i]    <- stats::quantile(x, probs[1])
    out$Upper[i]    <- stats::quantile(x, probs[2])
    out$exp_Estimate[i] <- exp(out$Estimate[i])
    out$exp_Lower[i]    <- exp(out$Lower[i])
    out$exp_Upper[i]    <- exp(out$Upper[i])
  }

  rownames(out) <- NULL
  out
}

#' Table of secondary attack rates
#'
#' Computes observed secondary attack rates (SAR) from the data, optionally
#' stratified by a covariate. Confidence intervals use the Wilson score method.
#'
#' @param fit An object of class \code{hhdynamics_fit}.
#' @param by Formula or character string naming the stratification variable
#'   (e.g. \code{~age} or \code{"age"}). Default: no stratification.
#' @return A data frame with columns: Stratum, N_contacts, N_infected, SAR,
#'   Lower, Upper.
#' @examples
#' \donttest{
#' data(inputdata)
#' fit <- household_dynamics(inputdata, ~sex, ~age,
#'   n_iteration = 15000, burnin = 5000, thinning = 1)
#' table_attack_rates(fit)
#' table_attack_rates(fit, by = ~age)
#' }
#' @export
table_attack_rates <- function(fit, by = NULL) {
  if (is.null(fit$input_data)) {
    stop("fit$input_data is NULL. Refit the model with household_dynamics() to use this function.",
         call. = FALSE)
  }

  dat <- fit$input_data
  contacts <- dat[dat$member > 0, ]

  if (is.null(by)) {
    # Overall SAR
    n_contacts <- nrow(contacts)
    n_infected <- sum(contacts$inf == 1)
    sar <- n_infected / n_contacts
    ci <- .wilson_ci(n_infected, n_contacts)
    return(data.frame(
      Stratum = "Overall",
      N_contacts = n_contacts,
      N_infected = n_infected,
      SAR = sar,
      Lower = ci[1],
      Upper = ci[2],
      stringsAsFactors = FALSE
    ))
  }

  # Parse by variable
  if (inherits(by, "formula")) {
    by_var <- all.vars(by)
  } else {
    by_var <- by
  }
  if (length(by_var) != 1) {
    stop("'by' must specify a single variable.", call. = FALSE)
  }
  if (!(by_var %in% names(contacts))) {
    stop(sprintf("Variable '%s' not found in the data.", by_var), call. = FALSE)
  }

  strata <- sort(unique(contacts[[by_var]]))
  strata <- strata[!is.na(strata)]  # drop NA strata
  n_strata <- length(strata)
  out <- data.frame(
    Stratum = as.character(strata),
    N_contacts = integer(n_strata),
    N_infected = integer(n_strata),
    SAR = numeric(n_strata),
    Lower = numeric(n_strata),
    Upper = numeric(n_strata),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_strata)) {
    sub <- contacts[which(contacts[[by_var]] == strata[i]), ]
    n_c <- nrow(sub)
    n_i <- sum(sub$inf == 1)
    ci <- .wilson_ci(n_i, n_c)
    out$N_contacts[i] <- n_c
    out$N_infected[i] <- n_i
    out$SAR[i] <- if (n_c > 0) n_i / n_c else NA
    out$Lower[i] <- ci[1]
    out$Upper[i] <- ci[2]
  }

  rownames(out) <- NULL
  out
}
