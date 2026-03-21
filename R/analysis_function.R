#' Run the MCMC for the household transmission model
#'
#' Low-level function that runs the C++ MCMC sampler. Most users should use
#' \code{\link{household_dynamics}} instead, which handles data preparation,
#' validation, and returns an S3 object with named parameters.
#'
#' @details
#' The MCMC uses a Metropolis-Hastings algorithm with adaptive proposal variances.
#' After 500 iterations, proposal standard deviations are set to the empirical
#' posterior standard deviation, with multiplicative tuning based on acceptance rate
#' (target: 20--30\%). Infection times for household contacts are jointly updated
#' at each iteration via a data augmentation step.
#'
#' The parameter vector has the following structure:
#' \enumerate{
#'   \item Standard deviation of random effect on infectivity (fixed at initial value if \code{with_rm = 0})
#'   \item Rate of infection from community (log scale)
#'   \item Rate of person-to-person transmission in households (log scale)
#'   \item Household size parameter (currently fixed at 0)
#'   \item Infectivity covariate effects (\code{n_inf} parameters)
#'   \item Susceptibility covariate effects (\code{n_sus} parameters)
#' }
#'
#' @param data_w The input data, in wide format (each row is a household), as produced by \code{\link{create_wide_data}}.
#' @param SI The mass function of the serial interval distribution. Defaults to the bundled influenza serial interval from Tsang et al. (2014).
#' @param n_iteration The number of iterations of the MCMC.
#' @param burnin The number of burn-in iterations to discard.
#' @param thinning The thinning interval for posterior samples.
#' @param n_inf The number of parameters affecting infectivity in the model.
#' @param n_sus The number of parameters affecting susceptibility in the model.
#' @param with_rm Indicator if the model has a random effect on individual infectivity (1) or not (0).
#' @param factor_group Integer vector mapping each dummy covariate column to its original factor group (from \code{\link{create_wide_data}}).
#' @param n_levels_vec Integer vector of the number of levels for each dummy column's factor (from \code{\link{create_wide_data}}).
#' @param estimate_SI Logical. If \code{TRUE}, jointly estimate Weibull shape/scale for the serial interval. Default is \code{FALSE}.
#' @return A list with 6 elements from the C++ MCMC:
#' \enumerate{
#'   \item Posterior samples matrix (post-burnin, thinned)
#'   \item Log-likelihood matrix (full chain, 3 columns: total, component 1, component 2)
#'   \item Random effect samples (post-burnin; empty if \code{with_rm = 0})
#'   \item Acceptance rates (per-parameter, numeric vector)
#'   \item Infection-time update acceptance rates (iterations x household members)
#'   \item Final imputed data matrix
#' }
#' @seealso \code{\link{household_dynamics}} for the high-level interface,
#'   \code{\link{create_wide_data}} for data preparation.
#' @examples
#' \donttest{
#' result_list <- create_wide_data(inputdata, ~sex, ~age)
#' data_w <- result_list[[1]]
#' n_inf <- result_list[[2]]
#' n_sus <- result_list[[3]]
#' mcmc_result <- run_MCMC(data_w,
#'   n_iteration = 15000, burnin = 5000,
#'   thinning = 1, n_inf = n_inf, n_sus = n_sus, with_rm = 0)
#' }
#' @export
run_MCMC <- function(data_w,SI = NULL,n_iteration = 15000,burnin = 5000,thinning = 1,n_inf,n_sus,with_rm,factor_group=integer(0),n_levels_vec=integer(0),estimate_SI=FALSE){
  if (is.null(SI)) {
    SI <- hhdynamics::SI
  }
  if (with_rm == 1) {
    warning("Random effects (with_rm = 1) is experimental. The random-effects output records one value per household (index case only), not per individual. Use with caution.", call. = FALSE)
  }

keep_iteration <- burnin + 1:((n_iteration - burnin)/thinning)*thinning
#### put to the MCMC
## first need to create parameter vector
# para 1. random effect on infectivity
# para 2. from community
# para 3. from household
# para 4. the household size parameter
# para about the input variables

para <-  c(1,0.01,0.1,0,rep(0.1,n_inf),rep(0.1,n_sus))

# append Weibull shape/scale if estimating SI
if (estimate_SI) {
  para <- c(para, 2.0, 4.0) # initial shape=2, scale=4
}

sigma <- c(1,rep(0.1,length(para)-1))
move <- rep(1,length(para))
move[1] <- with_rm
move[4] <- 0
sep1 <- 5
sep2 <- n_inf+n_sus+3 # 3 for onset time and the inf status, and the random effect

start_time <- Sys.time()
tt <- mcmc(as.matrix(data_w),SI,n_iteration,burnin,thinning,para,move,sigma,n_inf,n_sus,with_rm,sep1,sep2,as.integer(factor_group),as.integer(n_levels_vec),as.integer(estimate_SI))
end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time,start_time,units = "secs"))
message(paste0('The running time is ',round(elapsed), ' seconds'))

tt[[1]] <- tt[[1]][keep_iteration,]
attr(tt, "elapsed_time") <- elapsed
return(tt)
}

#################################################
## function to compute the mcmc output
## internal use, not output
para_summary <- function(mcmc_samples){
  y <- data.frame(
    mean = numeric(ncol(mcmc_samples)),
    lower = numeric(ncol(mcmc_samples)),
    upper = numeric(ncol(mcmc_samples)),
    acceptance = numeric(ncol(mcmc_samples))
  )
  for (i in 1:ncol(mcmc_samples)){
    q <- quantile(mcmc_samples[,i],c(0.025,0.975),na.rm=TRUE)
    y$mean[i] <- mean(mcmc_samples[,i],na.rm=TRUE)
    y$lower[i] <- q[1]
    y$upper[i] <- q[2]
    y$acceptance[i] <- sum(diff(mcmc_samples[,i])!=0)/nrow(mcmc_samples)
  }
  return(y)
}

#################################################
#' Create a wide format of household data
#'
#' Reshapes long-format household data (one row per individual) into wide format
#' (one row per household) for the C++ MCMC backend. Also constructs design
#' matrices for infectivity and susceptibility covariates.
#'
#' @details
#' The wide-format output has the following column structure per household:
#' \itemize{
#'   \item Columns 1--5: hhID, size, size (duplicate), index onset, end date
#'   \item Columns 6+: per-member blocks of (inf, onset, random effect placeholder,
#'     infectivity covariates, susceptibility covariates)
#' }
#' Missing values (members not present in a household) are filled with -1,
#' which the C++ code treats as a sentinel. Missing factor covariates are
#' coded as -99 (a distinct sentinel) and will be imputed during MCMC.
#'
#' @param input The input data, in long format. Must contain columns: hhID, member, inf, onset, size, end.
#' @param inf_factor Formula for factors affecting infectivity (e.g. \code{~sex} or \code{~sex + age}). Use \code{NULL} for no factors.
#' @param sus_factor Formula for factors affecting susceptibility (e.g. \code{~age}). Use \code{NULL} for no factors.
#' @return A list with 5 elements: (1) data frame in wide format, (2) number of infectivity parameters, (3) number of susceptibility parameters, (4) integer vector mapping each dummy column to its factor group, (5) integer vector of factor levels per dummy column.
#' @seealso \code{\link{household_dynamics}} for the high-level interface,
#'   \code{\link{run_MCMC}} for the low-level MCMC function.
#' @examples
#' wide_data <- create_wide_data(inputdata, ~sex, ~age)
#' @export
create_wide_data <- function(input,inf_factor,sus_factor){

  # first keep the essential
  data_a <- input[,c("hhID","member","inf","onset")]
  n_sus <- 0
  n_inf <- 0
  factor_group <- integer(0)
  n_levels_vec <- integer(0)
  group_counter <- 0L

  # get the infectivity
  data_a <- cbind(data_a,0)
  if (!is.null(inf_factor)){
    mf_inf <- stats::model.frame(inf_factor,input,na.action=stats::na.pass)
    design_matrix_inf <- model.matrix(attr(mf_inf,"terms"),mf_inf)
    n_inf <- ncol(design_matrix_inf)-1
    dm_inf <- design_matrix_inf[,-1,drop=FALSE]
    # Build factor group map for inf covariates
    assign_inf <- attr(design_matrix_inf,"assign")[-1]
    unique_terms_inf <- unique(assign_inf)
    for (ut in unique_terms_inf) {
      group_counter <- group_counter + 1L
      n_dummies <- sum(assign_inf == ut)
      nlev <- n_dummies + 1L
      factor_group <- c(factor_group, rep(group_counter, n_dummies))
      n_levels_vec <- c(n_levels_vec, rep(nlev, n_dummies))
    }
    # Replace NA with -99 sentinel (distinct from -1 padding)
    dm_inf[is.na(dm_inf)] <- -99
    data_a <- cbind(data_a,dm_inf)
  }
  # get the susceptibility
  if (!is.null(sus_factor)){
    mf_sus <- stats::model.frame(sus_factor,input,na.action=stats::na.pass)
    design_matrix_sus <- model.matrix(attr(mf_sus,"terms"),mf_sus)
    n_sus <- ncol(design_matrix_sus)-1
    dm_sus <- design_matrix_sus[,-1,drop=FALSE]
    # Build factor group map for sus covariates
    assign_sus <- attr(design_matrix_sus,"assign")[-1]
    unique_terms_sus <- unique(assign_sus)
    for (ut in unique_terms_sus) {
      group_counter <- group_counter + 1L
      n_dummies <- sum(assign_sus == ut)
      nlev <- n_dummies + 1L
      factor_group <- c(factor_group, rep(group_counter, n_dummies))
      n_levels_vec <- c(n_levels_vec, rep(nlev, n_dummies))
    }
    dm_sus[is.na(dm_sus)] <- -99
    data_a <- cbind(data_a,dm_sus)
  }


  data_w <- reshape(data_a,direction='wide',idvar='hhID',timevar = 'member')

  data_w <- merge(input[input$member==0,c("hhID","size","size","onset",'end')],data_w,by='hhID')

  # Set padding NAs to -1 (-99 sentinels for missing covariates are preserved)
  data_w[is.na(data_w)] <- -1

  return(list(data_w,n_inf,n_sus,as.integer(factor_group),as.integer(n_levels_vec)))
}


#################################################
#' Fit a household transmission model via MCMC
#'
#' The main function to fit the household transmission model to data. Estimates
#' the daily probability of infection from the community, the probability of
#' person-to-person transmission within households, and effects of covariates
#' on infectivity and susceptibility.
#'
#' @details
#' The model assumes that each household contact can be infected either from the
#' community (at a constant daily rate) or from an infected household member
#' (with probability governed by the serial interval distribution). Tertiary
#' transmission within households is accounted for. Infection times for
#' non-index cases are treated as latent variables and imputed via data
#' augmentation during MCMC.
#'
#' Covariate effects on infectivity and susceptibility enter multiplicatively
#' on the log scale. The \code{summary()} method reports exponentiated
#' estimates for interpretation as relative risks.
#'
#' \strong{Missing covariate imputation:} Factor covariates with missing values
#' (\code{NA}) are automatically imputed during MCMC via Bayesian data
#' augmentation, using a uniform categorical prior over factor levels. Only
#' factor covariates are supported; continuous covariates with \code{NA} will
#' produce an error. Interaction terms with missing data are not supported.
#'
#' The returned \code{hhdynamics_fit} object stores the full MCMC output,
#' enabling custom convergence diagnostics and post-processing. Key fields:
#' \describe{
#'   \item{\code{$samples}}{Posterior parameter samples (post-burnin, thinned). Columns named by parameter.}
#'   \item{\code{$log_likelihood}}{Log-likelihood trace for convergence assessment (full chain).}
#'   \item{\code{$acceptance}}{Per-parameter acceptance rates from the Metropolis-Hastings sampler.}
#'   \item{\code{$imputed_data}}{Final imputed dataset (wide format) with augmented infection times.}
#' }
#'
#' @param input The input data, in long format (each row is an individual). Required columns:
#'   \describe{
#'     \item{hhID}{Household identifier.}
#'     \item{member}{Member index (0 = index case, 1+ = contacts).}
#'     \item{size}{Number of individuals in the household.}
#'     \item{end}{End date of follow-up for that individual.}
#'     \item{inf}{Infection status (1 = infected, 0 = not). Index cases must have \code{inf = 1}.}
#'     \item{onset}{Onset time of symptoms.}
#'   }
#' @param inf_factor Formula for factors affecting infectivity (e.g. \code{~sex} or \code{~sex + age}). Use \code{NULL} for no factors. Default is \code{NULL}.
#' @param sus_factor Formula for factors affecting susceptibility (e.g. \code{~age}). Use \code{NULL} for no factors. Default is \code{NULL}.
#' @param SI The mass function of the serial interval distribution. Must be a numeric vector of length 14 summing to approximately 1. Defaults to the bundled influenza serial interval from Tsang et al. (2014). Not used when \code{estimate_SI = TRUE} (SI is estimated from data via Weibull parameterization).
#' @param n_iteration Total number of MCMC iterations. Default is 15000.
#' @param burnin Number of initial iterations to discard. Default is 5000.
#' @param thinning Thinning interval for posterior samples. Default is 1.
#' @param estimate_SI Logical. If \code{TRUE}, jointly estimate the serial interval distribution as a Weibull(shape, scale) alongside other model parameters. Two additional parameters (\code{si_shape}, \code{si_scale}) are added to the MCMC. Priors: shape ~ Uniform(0.1, 10), scale ~ Uniform(0.1, 20). Default is \code{FALSE}.
#' @return An object of class \code{\link{print.hhdynamics_fit}{hhdynamics_fit}}. Use \code{summary()} to get parameter estimates, \code{print()} for a brief overview, and \code{coef()} for posterior means. When \code{estimate_SI = TRUE}, the output includes \code{si_shape} and \code{si_scale} parameters.
#' @seealso \code{\link{simulate_data}} for simulating from the model,
#'   \code{\link{create_wide_data}} for data preparation,
#'   \code{\link{run_MCMC}} for the low-level MCMC interface.
#' @examples
#' \donttest{
#' data(inputdata)
#'
#' # Fit with default flu SI
#' fit <- household_dynamics(inputdata, ~sex, ~age,
#'   n_iteration = 15000, burnin = 5000, thinning = 1)
#' print(fit)
#' summary(fit)
#' coef(fit)
#'
#' # Fit without covariates (uses default SI)
#' fit2 <- household_dynamics(inputdata)
#' summary(fit2)
#'
#' # Jointly estimate SI from data
#' fit3 <- household_dynamics(inputdata, ~sex, ~age, estimate_SI = TRUE)
#' summary(fit3)  # includes si_shape and si_scale
#'
#' # Access MCMC samples for custom diagnostics
#' plot(fit$samples[, "community"], type = "l")
#' }
#' @export
household_dynamics <- function(input,inf_factor = NULL,sus_factor = NULL,SI = NULL,n_iteration = 15000,burnin = 5000,thinning = 1,estimate_SI = FALSE){
  if (is.null(SI)) {
    SI <- hhdynamics::SI
  }
  with_rm <- 0
  validate_inputs(input, inf_factor, sus_factor, SI, n_iteration, burnin, thinning, with_rm, estimate_SI)

  result_list <- create_wide_data(input,inf_factor,sus_factor)

  data_w <- result_list[[1]]
  n_inf <- result_list[[2]]
  n_sus <- result_list[[3]]
  factor_group <- result_list[[4]]
  n_levels_vec <- result_list[[5]]

  mcmc_result <- run_MCMC(data_w,SI,n_iteration,burnin,thinning,n_inf,n_sus,with_rm,factor_group,n_levels_vec,estimate_SI)

  # Build parameter names
  base_names <- c("re_sd", "community", "household", "size_param")
  if (n_inf + n_sus > 0){
    covariate_names <- names(data_w)[5 + 3 + (1:(n_inf + n_sus))]
  } else {
    covariate_names <- character(0)
  }
  param_names <- c(base_names, covariate_names)
  if (estimate_SI) {
    param_names <- c(param_names, "si_shape", "si_scale")
  }
  colnames(mcmc_result[[1]]) <- param_names

  # Determine which params get which transform in summary
  # 1=re_sd: none, 2-3=community/household: 1-exp(-x), 4=size: none, 5+=covariates: exp()
  param_transform <- rep("none", length(param_names))
  param_transform[2:3] <- "prob"
  n_base_plus_cov <- 4 + length(covariate_names)
  if (length(covariate_names) > 0) {
    param_transform[5:n_base_plus_cov] <- "exp"
  }
  # SI params stay "none" (already set)

  elapsed <- attr(mcmc_result, "elapsed_time")

  fit <- .new_hhdynamics_fit(
    samples = mcmc_result[[1]],
    log_likelihood = mcmc_result[[2]],
    random_effects = mcmc_result[[3]],
    acceptance = mcmc_result[[4]],
    update_accept = mcmc_result[[5]],
    imputed_data = mcmc_result[[6]],
    param_names = param_names,
    param_transform = param_transform,
    n_inf = n_inf,
    n_sus = n_sus,
    with_rm = with_rm,
    formula_inf = inf_factor,
    formula_sus = sus_factor,
    SI = SI,
    n_iteration = n_iteration,
    burnin = burnin,
    thinning = thinning,
    elapsed_time = elapsed,
    n_households = nrow(data_w),
    n_individuals = nrow(input)
  )

  return(fit)
}

#################################################
#' Simulate household transmission data
#'
#' Generates synthetic datasets from the household transmission model for
#' validation, power analysis, or posterior predictive checks.
#'
#' @details
#' The simulation uses the same household structure (sizes, follow-up periods,
#' covariate values) as the input data. The \code{rep_num} parameter replicates
#' the household structure to increase sample size. Infection outcomes and onset
#' times are simulated from the model given the parameter vector.
#'
#' The output is in wide format (one row per household), matching the internal
#' representation used by the C++ backend. Use this with \code{run_MCMC()} or
#' \code{household_dynamics()} to verify model recovery.
#'
#' @param input The dataset in long format (same structure as for \code{\link{household_dynamics}}).
#' @param rep_num The number of replications of the input dataset, to increase the sample size.
#' @param inf_factor Formula for factors affecting infectivity (e.g. \code{~sex}). Use \code{NULL} for no factors. Default is \code{NULL}.
#' @param sus_factor Formula for factors affecting susceptibility (e.g. \code{~age}). Use \code{NULL} for no factors. Default is \code{NULL}.
#' @param para The parameter vector, matching the structure from \code{\link{coef.hhdynamics_fit}}: (1) random effect SD, (2) community rate, (3) household rate, (4) size parameter, (5+) covariate effects.
#' @param SI The mass function of the serial interval distribution. Defaults to the bundled influenza serial interval from Tsang et al. (2014).
#' @param with_rm Indicator if the model has a random effect on individual infectivity (1) or not (0).
#' @return A simulated dataset in wide format (one row per household) based on the input parameter vectors.
#' @seealso \code{\link{household_dynamics}} for fitting the model,
#'   \code{\link{coef.hhdynamics_fit}} for extracting parameter estimates to use as simulation inputs.
#' @examples
#' \donttest{
#' data(inputdata)
#' data(SI)
#' para <- c(1, 0.01, 0.1, 0, 0.1, 0.1, 0.1)
#' simulated <- simulate_data(inputdata, 10, ~sex, ~age,
#'   SI = SI, para = para, with_rm = 0)
#' }
#' @export
simulate_data <- function(input,rep_num,inf_factor = NULL,sus_factor = NULL,SI = NULL,para,with_rm){
  if (is.null(SI)) {
    SI <- hhdynamics::SI
  }

  create_sim_data <- create_wide_data(input,inf_factor,sus_factor)
  simdata <- create_sim_data[[1]]

  n_inf <- create_sim_data[[2]]
  n_sus <- create_sim_data[[3]]

  # Force single-thread: R::runif inside parallelFor is not thread-safe,
  # causing correlated random draws and attenuated covariate effects
  prev_threads <- RcppParallel::defaultNumThreads()
  RcppParallel::setThreadOptions(numThreads = 1L)
  on.exit(RcppParallel::setThreadOptions(numThreads = prev_threads), add = TRUE)

  c1 <- sim_data(as.matrix(simdata[rep(1:nrow(simdata),rep_num),]),SI,para,n_inf,n_sus,with_rm,5,n_inf+n_sus+3)

  return(c1[[1]])
}
