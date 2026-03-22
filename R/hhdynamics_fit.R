# S3 class for hhdynamics model fit objects

#' @keywords internal
.new_hhdynamics_fit <- function(samples, log_likelihood, acceptance,
                                 update_accept, imputed_data, random_effects,
                                 param_names, param_transform, n_inf, n_sus,
                                 with_rm, formula_inf, formula_sus, SI,
                                 n_iteration, burnin, thinning, elapsed_time,
                                 n_households, n_individuals, input_data = NULL) {
  structure(list(
    samples = samples,
    log_likelihood = log_likelihood,
    acceptance = acceptance,
    update_accept = update_accept,
    imputed_data = imputed_data,
    random_effects = random_effects,
    param_names = param_names,
    param_transform = param_transform,
    n_inf = n_inf,
    n_sus = n_sus,
    with_rm = with_rm,
    formula_inf = formula_inf,
    formula_sus = formula_sus,
    SI = SI,
    n_iteration = n_iteration,
    burnin = burnin,
    thinning = thinning,
    elapsed_time = elapsed_time,
    n_households = n_households,
    n_individuals = n_individuals,
    input_data = input_data
  ), class = "hhdynamics_fit")
}

#' Print method for hhdynamics_fit
#'
#' @param x An object of class \code{hhdynamics_fit}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.hhdynamics_fit <- function(x, ...) {
  n_samples <- nrow(x$samples)
  cat("Household transmission model fit\n")
  cat(sprintf("  Data: %d households, %d individuals\n", x$n_households, x$n_individuals))
  cat(sprintf("  MCMC: %d post-burnin samples (%d iterations, burnin: %d, thin: %d)\n",
              n_samples, x$n_iteration, x$burnin, x$thinning))
  cat(sprintf("  Runtime: %s seconds\n", round(x$elapsed_time)))

  # Show which parameters are estimated
  show_params <- x$param_names
  show_params <- show_params[!show_params %in% c("re_sd", "size_param")]
  cat(sprintf("  Parameters: %s\n", paste(show_params, collapse = ", ")))

  cat("\nUse summary() for estimates, coef() for posterior means.\n")
  invisible(x)
}

#' Summary method for hhdynamics_fit
#'
#' Computes posterior summaries (mean, 2.5\%, 97.5\% credible intervals) for all model parameters.
#' For community and household parameters, reports the daily probability (via 1-exp(-x) transform).
#' For covariate effects, additionally reports exponentiated estimates.
#'
#' @param object An object of class \code{hhdynamics_fit}.
#' @param ... Additional arguments (unused).
#' @return An object of class \code{summary.hhdynamics_fit}, which is a data frame with columns:
#'   Variable, Point.estimate, Lower.bound, Upper.bound, and optionally exp columns for covariates.
#' @export
summary.hhdynamics_fit <- function(object, ...) {
  ps <- para_summary(object$samples)

  output <- data.frame(matrix(NA, nrow(ps), 7))
  names(output) <- c("Variable", "Point estimate", "Lower bound", "Upper bound",
                      "exp(Point estimate)", "exp(Lower bound)", "exp(Upper bound)")

  z1 <- as.matrix(ps[, c("mean", "lower", "upper", "acceptance")])

  # Transform community and household params to probability scale
  z1[2:3, 1:3] <- 1 - exp(-z1[2:3, 1:3])

  z1_exp <- exp(z1[, 1:3, drop = FALSE])
  z1_combined <- cbind(z1, z1_exp)

  output[] <- z1_combined[, c(4, 1:3, 5:7)]

  # Variable names
  var_names <- c("Standard deviation of random effect",
                 "Daily probability of infection from community",
                 "Probability of person-to-person transmission in households",
                 "Parameter of the relationship between transmission and number of household members")
  if (length(object$param_names) > 4) {
    var_names <- c(var_names, object$param_names[5:length(object$param_names)])
  }
  output[, 1] <- var_names

  # Clear exp columns for base params (not meaningful)
  output[1:4, 5:7] <- NA

  # Clear exp columns for SI params (reported on natural scale)
  si_rows <- which(object$param_names %in% c("si_shape", "si_scale"))
  if (length(si_rows) > 0) {
    output[si_rows, 5:7] <- NA
  }

  # Drop size parameter row (always fixed at 0)
  output <- output[-4, ]

  # Drop RE row if no random effect
  if (object$with_rm == 0) {
    output <- output[-1, ]
  }

  rownames(output) <- NULL

  structure(output, class = c("summary.hhdynamics_fit", "data.frame"))
}

#' Print method for summary.hhdynamics_fit
#'
#' @param x An object of class \code{summary.hhdynamics_fit}.
#' @param digits Number of significant digits for printing. Default is 3.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.summary.hhdynamics_fit <- function(x, digits = 3, ...) {
  output_print <- x
  output_print[, -1] <- round(output_print[, -1], digits)
  print.data.frame(output_print, row.names = FALSE)
  invisible(x)
}

#' Extract model coefficients from hhdynamics_fit
#'
#' Returns named vector of posterior means for all estimated parameters.
#'
#' @param object An object of class \code{hhdynamics_fit}.
#' @param ... Additional arguments (unused).
#' @return A named numeric vector of posterior means.
#' @export
coef.hhdynamics_fit <- function(object, ...) {
  means <- colMeans(object$samples)
  names(means) <- object$param_names
  means
}
