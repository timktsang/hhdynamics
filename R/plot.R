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

#' Plot secondary attack rates
#'
#' Forest plot showing observed secondary attack rates (SAR) with Wilson score
#' 95\% confidence intervals. Supports overall SAR, stratification by one or
#' more covariates, and combinations thereof. Variable names are used as bold
#' section headers; strata are indented below. The layout mirrors
#' \code{plot_covariates()}: labels on the left, point estimates and CI bars in
#' the middle, and \code{n/N} counts plus SAR percentages on the right.
#'
#' @param fit An object of class \code{hhdynamics_fit}.
#' @param by Formula, character string, or a \emph{list} of formulas naming the
#'   stratification variable(s). Examples: \code{~age}, \code{"age"},
#'   \code{list(~sex, ~age)}. Default: \code{NULL} (overall SAR only).
#' @param include_overall Logical. When \code{TRUE}, an "Overall" row is
#'   prepended even when \code{by} is specified. Default: \code{FALSE}.
#' @param labels Optional named list of custom display labels for variables
#'   and their levels. Each element is a list with \code{name} (display name
#'   for the section header) and/or \code{levels} (character vector of level
#'   labels in the same order as \code{sort(unique(variable))}). Names must
#'   match variable names in the data. Example:
#'   \code{list(age = list(name = "Age Group", levels = c("0-5", "6-17", "18+")))}.
#' @param file Optional file path for PDF output. Height is auto-calculated
#'   from the number of rows. Default: \code{NULL} (current device).
#' @param width PDF width in inches. Default: 8.
#' @param height PDF height in inches. Default: \code{0.45 * n_rows + 1.8}.
#' @param xlim Numeric vector of length 2 for the x-axis range (probability
#'   scale). Default: auto-determined from the data.
#' @param cex Character expansion factor. Default: 0.85.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#' @return Invisible data frame of the estimate rows (Stratum, N_contacts,
#'   N_infected, SAR, Lower, Upper).
#' @examples
#' \donttest{
#' data(inputdata)
#' fit <- household_dynamics(inputdata, ~sex, ~age,
#'   n_iteration = 15000, burnin = 5000, thinning = 1)
#'
#' # Overall only
#' plot_attack_rate(fit)
#'
#' # Stratified by age with section header
#' plot_attack_rate(fit, by = ~age)
#'
#' # Combined: overall + age + sex in one figure
#' plot_attack_rate(fit, by = list(~sex, ~age), include_overall = TRUE,
#'   labels = list(sex = list(name = "Sex", levels = c("Male", "Female")),
#'                 age = list(name = "Age Group", levels = c("0-5", "6-17", "18+"))))
#' }
#' @export
plot_attack_rate <- function(fit, by = NULL, include_overall = FALSE,
                             labels = NULL,
                             file = NULL, width = 8, height = NULL,
                             xlim = NULL, cex = 0.85, ...) {

  # --- Build row list -------------------------------------------------------
  rows <- list()

  # Overall row (always shown when by = NULL, or when include_overall = TRUE)
  if (is.null(by) || include_overall) {
    tab_ov <- table_attack_rates(fit, by = NULL)
    rows[[length(rows) + 1]] <- list(
      type = "estimate", label = "Overall", indented = FALSE,
      SAR        = tab_ov$SAR,
      Lower      = tab_ov$Lower,
      Upper      = tab_ov$Upper,
      N_infected = tab_ov$N_infected,
      N_contacts = tab_ov$N_contacts
    )
  }

  # Per-variable sections
  if (!is.null(by)) {
    if (inherits(by, "formula") || is.character(by)) by <- list(by)

    for (f in by) {
      var_name <- if (inherits(f, "formula")) all.vars(f) else f
      if (length(var_name) != 1)
        stop("Each 'by' entry must specify a single variable.", call. = FALSE)

      # Resolve display name and level labels from labels argument
      var_display   <- var_name
      custom_levels <- NULL
      if (!is.null(labels) && var_name %in% names(labels)) {
        lab <- labels[[var_name]]
        if (!is.null(lab$name))   var_display   <- lab$name
        if (!is.null(lab$levels)) custom_levels <- lab$levels
      }

      # Section header row
      rows[[length(rows) + 1]] <- list(type = "header", label = var_display)

      # Estimate rows
      tab <- table_attack_rates(fit, by = if (inherits(f, "formula")) f else stats::as.formula(paste0("~", f)))
      contacts   <- fit$input_data[fit$input_data$member > 0, ]
      data_levs  <- as.character(sort(unique(contacts[[var_name]][!is.na(contacts[[var_name]])])))

      for (i in seq_len(nrow(tab))) {
        strat_lbl <- as.character(tab$Stratum[i])
        if (!is.null(custom_levels)) {
          idx <- match(strat_lbl, data_levs)
          if (!is.na(idx) && idx <= length(custom_levels))
            strat_lbl <- custom_levels[idx]
        }
        rows[[length(rows) + 1]] <- list(
          type = "estimate", label = strat_lbl, indented = TRUE,
          SAR        = tab$SAR[i],
          Lower      = tab$Lower[i],
          Upper      = tab$Upper[i],
          N_infected = tab$N_infected[i],
          N_contacts = tab$N_contacts[i]
        )
      }
    }
  }

  n_rows <- length(rows)

  # --- Open PDF device -------------------------------------------------------
  opened_device <- FALSE
  if (!is.null(file)) {
    if (is.null(height)) height <- 0.45 * n_rows + 1.8
    grDevices::pdf(file, width = width, height = height)
    opened_device <- TRUE
  }

  # --- X-axis range ----------------------------------------------------------
  est_rows <- Filter(function(r) r$type == "estimate", rows)
  all_upper <- sapply(est_rows, function(r) r$Upper)
  if (is.null(xlim)) {
    x_max <- min(1, max(all_upper, na.rm = TRUE) * 1.15)
    x_max <- max(x_max, 0.2)
    xlim <- c(0, x_max)
  }
  dx <- diff(xlim)

  # Coordinate layout: [label region] | [forest 0..xlim[2]] | [n/N] [SAR text]
  label_x    <- xlim[1] - dx
  count_x    <- xlim[2] + 0.06 * dx
  right_edge <- xlim[2] + 0.70 * dx

  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(op)
    if (opened_device) grDevices::dev.off()
  })
  graphics::par(mar = c(0, 0, 0, 0))

  y_top    <- n_rows + 1.5
  y_bottom <- -1.8
  graphics::plot(0, 0, type = "n",
                 xlim = c(label_x, right_edge),
                 ylim = c(y_bottom, y_top),
                 axes = FALSE, xlab = "", ylab = "", ...)

  # Column headers
  header_y <- n_rows + 0.8
  graphics::text(label_x + 0.04 * dx, header_y, "Stratum",
                 adj = 0, font = 2, cex = cex)
  graphics::text(count_x, header_y, "n / N",
                 adj = 0, font = 2, cex = cex)
  graphics::text(right_edge, header_y, "SAR",
                 adj = 1, font = 2, cex = cex)
  graphics::text(right_edge, header_y - 0.55, "(95% CI)",
                 adj = 1, cex = cex * 0.85)

  # --- Draw rows -------------------------------------------------------------
  for (i in seq_len(n_rows)) {
    row <- rows[[i]]
    y <- n_rows - i + 0.5

    if (i %% 2 == 0) {
      graphics::rect(label_x, y - 0.5, right_edge, y + 0.5,
                     col = grDevices::rgb(0, 0, 0, 0.05), border = NA)
    }

    if (row$type == "header") {
      graphics::text(label_x + 0.04 * dx, y, row$label,
                     adj = 0, font = 2, cex = cex * 0.95)

    } else {
      indent_x <- if (isTRUE(row$indented)) label_x + 0.16 * dx else label_x + 0.04 * dx
      graphics::text(indent_x, y, row$label, adj = 0, cex = cex)

      if (!is.na(row$SAR)) {
        lo <- max(xlim[1], row$Lower)
        hi <- min(xlim[2], row$Upper)
        pt <- max(xlim[1], min(xlim[2], row$SAR))

        graphics::points(pt, y, pch = 15, cex = 0.9)
        graphics::segments(lo, y, hi, y, lwd = 1.5)

        if (row$Lower < xlim[1])
          graphics::arrows(lo + 0.04 * dx, y, lo, y, length = 0.05)
        if (row$Upper > xlim[2])
          graphics::arrows(hi - 0.04 * dx, y, hi, y, length = 0.05)
      }

      graphics::text(count_x, y,
                     sprintf("%d / %d", as.integer(row$N_infected), as.integer(row$N_contacts)),
                     adj = 0, cex = cex * 0.9)

      sar_text <- if (!is.na(row$SAR)) {
        sprintf("%.1f%% (%.1f, %.1f)",
                row$SAR * 100, row$Lower * 100, row$Upper * 100)
      } else "\u2014"
      graphics::text(right_edge, y, sar_text, adj = 1, cex = cex * 0.9)
    }
  }

  # Dashed reference line at x = 0 (spanning estimate rows only)
  est_idx <- which(sapply(rows, function(r) r$type == "estimate"))
  if (length(est_idx) > 0) {
    graphics::segments(0, n_rows - max(est_idx), 0, n_rows - min(est_idx) + 1,
                       lty = 2, col = "grey60")
  }

  # X-axis with percentage labels
  axis_y    <- -0.1
  tick_vals <- pretty(xlim, n = 5)
  tick_vals <- tick_vals[tick_vals >= xlim[1] & tick_vals <= xlim[2]]
  graphics::axis(1, at = tick_vals,
                 labels = paste0(round(tick_vals * 100), "%"),
                 pos = axis_y, cex.axis = cex)

  # Return estimate rows as a data frame
  out <- do.call(rbind, lapply(est_rows, function(r) {
    data.frame(Stratum = r$label, N_contacts = r$N_contacts,
               N_infected = r$N_infected, SAR = r$SAR,
               Lower = r$Lower, Upper = r$Upper, stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  invisible(out)
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
