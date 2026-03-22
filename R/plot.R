# Plot functions for hhdynamics

#' Plot method for hhdynamics_fit objects
#'
#' @param x An object of class \code{hhdynamics_fit}.
#' @param type Type of plot: \code{"diagnostics"} (default), \code{"transmission"},
#'   \code{"attack_rate"}, or \code{"covariates"}.
#' @param ... Additional arguments passed to the underlying plot function.
#' @return Invisible NULL.
#' @export
plot.hhdynamics_fit <- function(x, type = "diagnostics", ...) {
  type <- match.arg(type, c("diagnostics", "transmission", "attack_rate", "covariates"))
  switch(type,
    diagnostics  = plot_diagnostics(x, ...),
    transmission = plot_transmission(x, ...),
    attack_rate  = plot_attack_rate(x, ...),
    covariates   = plot_covariates(x, ...)
  )
}

#' MCMC diagnostic plots
#'
#' Produces trace plots and posterior density plots for each estimated parameter.
#' Trace plots show the MCMC chain with the posterior mean (red dashed line).
#' Density plots show the marginal posterior with 95\% credible interval bounds
#' (blue dashed lines) and effective sample size (ESS) in the title.
#'
#' Community and household parameters are shown on the probability scale
#' (via the \code{1 - exp(-x)} transform).
#'
#' @param fit An object of class \code{hhdynamics_fit}.
#' @param params Optional character vector of parameter names to plot. If
#'   \code{NULL} (default), all estimated parameters are plotted (fixed
#'   parameters like \code{size_param} are skipped).
#' @return Invisible NULL. Called for its side effect of producing plots.
#' @examples
#' \donttest{
#' data(inputdata)
#' fit <- household_dynamics(inputdata, ~sex, ~age,
#'   n_iteration = 15000, burnin = 5000, thinning = 1)
#' plot_diagnostics(fit)
#' plot_diagnostics(fit, params = c("community", "household"))
#' }
#' @export
plot_diagnostics <- function(fit, params = NULL) {
  plot_names <- .get_plot_params(fit, params)
  n_params <- length(plot_names)
  idx <- match(plot_names, fit$param_names)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mfrow = c(n_params, 2),
                mar = c(3, 3, 2.5, 1),
                mgp = c(2, 0.7, 0))

  for (i in seq_len(n_params)) {
    x_raw <- fit$samples[, idx[i]]
    transform <- fit$param_transform[idx[i]]
    pname <- plot_names[i]

    # Display label
    display_name <- switch(pname,
      "community" = "Community (probability)",
      "household" = "Household (probability)",
      pname
    )

    # Transform for density (community/household to probability scale)
    x_display <- .transform_param(x_raw, transform)

    # ESS computed on raw scale
    ess <- round(.effective_sample_size(x_raw))

    # Trace plot (raw scale)
    graphics::plot(x_raw, type = "l",
                   main = pname,
                   xlab = "Iteration (post-burnin)",
                   ylab = "Value",
                   col = "grey40")
    graphics::abline(h = mean(x_raw), col = "red", lty = 2)

    # Density plot (transformed scale)
    d <- stats::density(x_display)
    graphics::plot(d, main = paste0(display_name, "  (ESS: ", ess, ")"),
                   xlab = display_name,
                   ylab = "Density")
    ci <- stats::quantile(x_display, c(0.025, 0.975))
    graphics::abline(v = ci, col = "blue", lty = 2)
    graphics::abline(v = mean(x_display), col = "red", lty = 2)
  }

  invisible(NULL)
}

#' Plot transmission probability over time since onset
#'
#' Shows the daily probability of person-to-person transmission as a function
#' of days since the infector's symptom onset. The serial interval distribution
#' shapes this curve. The median and 95\% credible interval are computed from
#' the posterior samples.
#'
#' When the model was fitted with \code{estimate_SI = TRUE}, the uncertainty
#' band incorporates serial interval uncertainty (via the Weibull shape/scale
#' posterior).
#'
#' @param fit An object of class \code{hhdynamics_fit}.
#' @param hh_size Reference household size. Default: median from data.
#' @param col Color for the median line and credible interval polygon.
#'   Default: \code{"steelblue"}.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#' @return Invisible data frame with columns Day, Median, Lower, Upper.
#' @examples
#' \donttest{
#' data(inputdata)
#' fit <- household_dynamics(inputdata, ~sex, ~age,
#'   n_iteration = 15000, burnin = 5000, thinning = 1)
#' plot_transmission(fit)
#' }
#' @export
plot_transmission <- function(fit, hh_size = NULL, col = "steelblue", ...) {
  samples <- fit$samples
  n_draws <- nrow(samples)
  si_len <- length(fit$SI)
  days <- seq_len(si_len)

  # Reference household size
  if (is.null(hh_size)) {
    if (!is.null(fit$input_data)) {
      hh_size <- stats::median(fit$input_data$size[fit$input_data$member == 0])
    } else {
      hh_size <- 4  # fallback
    }
  }

  household_rate <- samples[, "household"]
  size_param <- samples[, "size_param"]

  # Compute SI for each draw
  has_si_est <- "si_shape" %in% colnames(samples)
  si_matrix <- matrix(NA, n_draws, si_len)
  if (has_si_est) {
    shapes <- samples[, "si_shape"]
    scales <- samples[, "si_scale"]
    for (t in seq_len(si_len)) {
      si_matrix[, t] <- stats::pweibull(t + 1, shapes, scales) -
                         stats::pweibull(t, shapes, scales)
    }
  } else {
    si_matrix <- matrix(fit$SI, nrow = n_draws, ncol = si_len, byrow = TRUE)
  }

  # P(infection on day t) = 1 - exp(-beta_hh * SI[t] / (hh_size - 1)^alpha)
  prob_matrix <- matrix(NA, n_draws, si_len)
  for (t in seq_len(si_len)) {
    hazard <- household_rate * si_matrix[, t] / (hh_size - 1)^size_param
    prob_matrix[, t] <- 1 - exp(-hazard)
  }

  # Summarize
  med <- apply(prob_matrix, 2, stats::median)
  lower <- apply(prob_matrix, 2, stats::quantile, 0.025)
  upper <- apply(prob_matrix, 2, stats::quantile, 0.975)

  result <- data.frame(Day = days, Median = med, Lower = lower, Upper = upper)

  # Plot
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mar = c(4, 4, 2, 1))

  y_max <- max(upper) * 1.1
  graphics::plot(days, med, type = "n",
                 xlim = c(1, si_len), ylim = c(0, y_max),
                 xlab = "Days since infector onset",
                 ylab = "P(transmission per day)",
                 main = paste0("Household transmission probability (size = ", hh_size, ")"),
                 ...)
  # CI polygon
  graphics::polygon(c(days, rev(days)), c(lower, rev(upper)),
                    col = grDevices::adjustcolor(col, alpha.f = 0.3), border = NA)
  graphics::lines(days, med, col = col, lwd = 2)
  graphics::points(days, med, col = col, pch = 16, cex = 0.8)

  invisible(result)
}

#' Plot secondary attack rates by covariate
#'
#' Bar plot showing observed secondary attack rates (SAR) with Wilson score
#' 95\% confidence intervals, optionally stratified by a covariate.
#'
#' @param fit An object of class \code{hhdynamics_fit}.
#' @param by Formula or character string naming the stratification variable
#'   (e.g. \code{~age} or \code{"age"}). Default: no stratification (overall SAR).
#' @param col Colors for bars. Default: grey palette.
#' @param ... Additional graphical parameters.
#' @return Invisible data frame with the plotted attack rate values.
#' @examples
#' \donttest{
#' data(inputdata)
#' fit <- household_dynamics(inputdata, ~sex, ~age,
#'   n_iteration = 15000, burnin = 5000, thinning = 1)
#' plot_attack_rate(fit)
#' plot_attack_rate(fit, by = ~age)
#' }
#' @export
plot_attack_rate <- function(fit, by = NULL, col = NULL, ...) {
  tab <- table_attack_rates(fit, by = by)

  if (is.null(col)) {
    col <- grDevices::grey.colors(nrow(tab), start = 0.4, end = 0.7)
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mar = c(5, 4, 3, 1))

  y_max <- min(1, max(tab$Upper, na.rm = TRUE) * 1.2)
  bp <- graphics::barplot(tab$SAR, names.arg = tab$Stratum,
                          ylim = c(0, y_max),
                          col = col,
                          ylab = "Secondary attack rate",
                          main = "Observed secondary attack rates",
                          ...)
  # Error bars
  graphics::arrows(bp, tab$Lower, bp, tab$Upper,
                   angle = 90, code = 3, length = 0.05)

  invisible(tab)
}

#' Plot household infection timeline
#'
#' Visualizes the infection timeline for a single household. The index case
#' is shown as a filled triangle, infected contacts as filled circles at
#' their (imputed) onset times, and uninfected contacts as open circles
#' spanning their follow-up period.
#'
#' @param fit An object of class \code{hhdynamics_fit}.
#' @param hh_id Household identifier to visualize.
#' @param col Colors for infected and uninfected members.
#'   Default: \code{c("firebrick", "grey60")}.
#' @param ... Additional graphical parameters.
#' @return Invisible NULL.
#' @examples
#' \donttest{
#' data(inputdata)
#' fit <- household_dynamics(inputdata, ~sex, ~age,
#'   n_iteration = 15000, burnin = 5000, thinning = 1)
#' plot_household(fit, hh_id = 1)
#' }
#' @export
plot_household <- function(fit, hh_id, col = NULL, ...) {
  if (is.null(fit$input_data)) {
    stop("fit$input_data is NULL. Refit the model with household_dynamics() to use this function.",
         call. = FALSE)
  }

  if (is.null(col)) col <- c("firebrick", "grey60")

  hh <- fit$input_data[fit$input_data$hhID == hh_id, ]
  if (nrow(hh) == 0) {
    stop(sprintf("Household ID %s not found in the data.", hh_id), call. = FALSE)
  }
  hh <- hh[order(hh$member), ]

  n_members <- nrow(hh)
  idx_onset <- hh$onset[hh$member == 0]
  end_time <- hh$end[1]

  # Get imputed onset times from the final MCMC iteration
  sep1 <- 5
  sep2 <- fit$n_inf + fit$n_sus + 3
  imp <- fit$imputed_data
  # Find matching household row
  hh_rows <- which(imp[, 1] == hh_id)
  if (length(hh_rows) == 0) {
    # Fallback: use input data onsets
    imputed_onsets <- hh$onset
  } else {
    hh_row <- hh_rows[1]
    imputed_onsets <- numeric(n_members)
    for (m in 0:(n_members - 1)) {
      imputed_onsets[m + 1] <- imp[hh_row, sep1 + m * sep2 + 2]  # R 1-indexed onset
    }
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mar = c(4, 6, 3, 1))

  x_range <- c(min(idx_onset, na.rm = TRUE) - 1, end_time + 1)
  graphics::plot(NULL, xlim = x_range, ylim = c(0.5, n_members + 0.5),
                 xlab = "Day", ylab = "",
                 yaxt = "n",
                 main = paste("Household", hh_id),
                 ...)
  graphics::axis(2, at = seq_len(n_members),
                 labels = paste("Member", hh$member), las = 1)

  for (i in seq_len(n_members)) {
    y <- i
    is_infected <- hh$inf[i] == 1
    is_index <- hh$member[i] == 0
    onset_val <- imputed_onsets[i]

    # Follow-up line
    graphics::segments(idx_onset, y, end_time, y, col = "grey80", lwd = 1)

    if (is_index) {
      graphics::points(onset_val, y, pch = 17, col = col[1], cex = 1.5)  # triangle
    } else if (is_infected) {
      graphics::points(onset_val, y, pch = 16, col = col[1], cex = 1.3)  # filled circle
    } else {
      graphics::points(end_time, y, pch = 1, col = col[2], cex = 1.0)  # open circle at end
    }
  }

  graphics::legend("topright",
                   legend = c("Index case", "Infected contact", "Uninfected"),
                   pch = c(17, 16, 1),
                   col = c(col[1], col[1], col[2]),
                   bty = "n", cex = 0.8)

  invisible(NULL)
}

#' Forest plot of covariate effects
#'
#' Produces a forest plot showing estimated relative risks for covariate effects
#' on susceptibility and infectiousness. Covariates are grouped by variable with
#' bold headers, reference categories labeled, alternating row shading, and
#' estimate text with credible intervals on the right.
#'
#' When \code{file} is provided, the plot is saved to a PDF with dimensions
#' automatically calculated from the number of covariate rows. When \code{file}
#' is \code{NULL}, the plot is drawn to the current graphics device.
#'
#' @param fit An object of class \code{hhdynamics_fit}.
#' @param probs Numeric vector of length 2 for credible interval bounds.
#'   Default: \code{c(0.025, 0.975)} (95\% CrI).
#' @param labels Optional named list of custom labels for covariates. Each
#'   element is a list with \code{name} (display name for the variable header)
#'   and \code{levels} (character vector of level labels, including the
#'   reference level first). Names must match variable names in the formula.
#'   Example: \code{list(sex = list(name = "Sex", levels = c("Male", "Female")))}.
#' @param file Optional file path for PDF output. When provided, a PDF is
#'   created with auto-calculated width and height based on the number of rows.
#'   Default: \code{NULL} (plot to current device).
#' @param width PDF width in inches. Default: 11. Only used when \code{file}
#'   is not \code{NULL}.
#' @param height PDF height in inches. Default: auto-calculated as
#'   \code{0.45 * n_rows + 1.8}. Only used when \code{file} is not \code{NULL}.
#' @param xlim Numeric vector of length 2 for the x-axis range on the natural
#'   scale. Default: auto-determined from the credible intervals.
#' @param xlab_left Label for the left direction arrow. Default: \code{"Lower Risk"}.
#' @param xlab_right Label for the right direction arrow. Default: \code{"Higher Risk"}.
#' @param cex Character expansion factor. Default: 0.85.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#' @return Invisible NULL. Called for its side effect of producing a plot.
#' @examples
#' \donttest{
#' data(inputdata)
#' fit <- household_dynamics(inputdata, ~sex, ~age,
#'   n_iteration = 15000, burnin = 5000, thinning = 1)
#' plot_covariates(fit)
#'
#' # Save to PDF with auto-sized dimensions
#' plot_covariates(fit, file = "covariates.pdf",
#'   labels = list(sex = list(name = "Sex", levels = c("Male", "Female")),
#'                 age = list(name = "Age Group", levels = c("0-5", "6-17", "18+"))))
#' }
#' @export
plot_covariates <- function(fit, probs = c(0.025, 0.975),
                            labels = NULL,
                            file = NULL,
                            width = 11,
                            height = NULL,
                            xlim = NULL,
                            xlab_left = "Lower Risk",
                            xlab_right = "Higher Risk",
                            cex = 0.85, ...) {
  n_cov <- fit$n_inf + fit$n_sus
  if (n_cov == 0) {
    message("No covariates in the model.")
    return(invisible(NULL))
  }
  if (is.null(fit$input_data)) {
    stop("fit$input_data is NULL. Refit with household_dynamics() to use this function.",
         call. = FALSE)
  }

  rows <- .build_forest_rows(fit, probs, labels)
  n_rows <- length(rows)

  # Auto-size and open PDF device if file is specified
  opened_device <- FALSE
  if (!is.null(file)) {
    if (is.null(height)) {
      height <- 0.45 * n_rows + 1.8
    }
    grDevices::pdf(file, width = width, height = height)
    opened_device <- TRUE
  }

  # Determine x-axis range (log2 scale) from estimate rows
  est_rows <- Filter(function(r) r$type == "estimate", rows)
  all_vals <- unlist(lapply(est_rows, function(r) c(r$lower, r$upper)))
  if (is.null(xlim)) {
    log_min <- floor(log2(min(all_vals, na.rm = TRUE)))
    log_max <- ceiling(log2(max(all_vals, na.rm = TRUE)))
    log_min <- max(log_min, -4)
    log_max <- min(log_max, 4)
    xlim <- 2^c(log_min, log_max)
  }
  log_xlim <- log2(xlim)

  # Coordinate layout:
  # label region | forest plot (log2 scale) | estimate text
  label_x <- log_xlim[1] - 8
  text_x <- log_xlim[2] + 0.8
  right_edge <- text_x + 4.5

  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(op)
    if (opened_device) grDevices::dev.off()
  })
  graphics::par(mar = c(0, 0, 0, 0))

  # Compact vertical layout: 1 unit per row
  y_top <- n_rows + 1.5
  y_bottom <- -2.0
  graphics::plot(0, 0, type = "n",
                 xlim = c(label_x, right_edge),
                 ylim = c(y_bottom, y_top),
                 axes = FALSE, xlab = "", ylab = "", ...)

  # Column headers (single row)
  header_y <- n_rows + 0.8
  graphics::text(label_x + 0.3, header_y, "Characteristics",
                 adj = 0, font = 2, cex = cex)
  graphics::text(right_edge - 0.3, header_y, "Risk Ratio",
                 adj = 1, font = 2, cex = cex)
  graphics::text(right_edge - 0.3, header_y - 0.6, "(95% CrI)",
                 adj = 1, cex = cex * 0.85)

  # Draw each row (top to bottom, compact spacing)
  for (i in seq_len(n_rows)) {
    row <- rows[[i]]
    y <- n_rows - i + 0.5

    # Alternating row shading
    if (i %% 2 == 0) {
      graphics::rect(label_x, y - 0.5, right_edge, y + 0.5,
                     col = grDevices::rgb(0, 0, 0, 0.05), border = NA)
    }

    if (row$type == "section") {
      graphics::text(label_x + 0.3, y, row$label,
                     adj = 0, font = 2, cex = cex)

    } else if (row$type == "header") {
      graphics::text(label_x + 1, y, row$label,
                     adj = 0, font = 2, cex = cex * 0.95)

    } else if (row$type == "ref") {
      graphics::text(label_x + 2, y, row$label,
                     adj = 0, font = 1, cex = cex)
      graphics::text(right_edge - 0.3, y, "Ref",
                     adj = 1, cex = cex)

    } else if (row$type == "estimate") {
      # Label — factor levels indented more than continuous covariates
      indent_x <- if (row$is_level) label_x + 2 else label_x + 1
      graphics::text(indent_x, y, row$label,
                     adj = 0, font = 1, cex = cex)

      # Point estimate and CI line
      log_est <- log2(row$estimate)
      log_lo <- log2(max(row$lower, xlim[1]))
      log_hi <- log2(min(row$upper, xlim[2]))

      graphics::points(log_est, y, pch = 16, cex = 0.8)
      graphics::segments(log_lo, y, log_hi, y)

      # Arrows if CI extends beyond plot limits
      if (row$lower < xlim[1]) {
        graphics::arrows(log_lo + 0.3, y, log_lo, y, length = 0.05)
      }
      if (row$upper > xlim[2]) {
        graphics::arrows(log_hi - 0.3, y, log_hi, y, length = 0.05)
      }

      # Estimate text (right-aligned)
      est_text <- sprintf("%.2f (%.2f, %.2f)",
                          row$estimate, row$lower, row$upper)
      graphics::text(right_edge - 0.3, y, est_text,
                     adj = 1, cex = cex * 0.9)
    }
  }

  # Dashed vertical reference line at 1 (log2(1) = 0)
  y_est_idx <- which(sapply(rows, function(r) r$type == "estimate"))
  if (length(y_est_idx) > 0) {
    y_line_top <- n_rows - min(y_est_idx) + 1
    y_line_bot <- n_rows - max(y_est_idx)
    graphics::segments(0, y_line_bot, 0, y_line_top, lty = 2)
  }

  # X-axis (log2 scale with natural-scale labels)
  axis_y <- -0.1
  ticks <- seq(log_xlim[1], log_xlim[2])
  graphics::axis(1, at = ticks, labels = 2^ticks, pos = axis_y, cex.axis = cex)

  # Direction arrows and labels (below axis with clearance)
  arrow_y <- axis_y - 1.0
  graphics::arrows(-0.15, arrow_y, log_xlim[1], arrow_y, length = 0.05)
  graphics::arrows(0.15, arrow_y, log_xlim[2], arrow_y, length = 0.05)
  graphics::text(log_xlim[1] / 2, arrow_y - 0.3, xlab_left, cex = cex * 0.85)
  graphics::text(log_xlim[2] / 2, arrow_y - 0.3, xlab_right, cex = cex * 0.85)

  invisible(NULL)
}
