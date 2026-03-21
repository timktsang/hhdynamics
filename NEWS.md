# hhdynamics 1.2.0

## New features

* **Optional serial interval estimation.** Set `estimate_SI = TRUE` in
  `household_dynamics()` to jointly estimate the serial interval distribution
  as a Weibull(shape, scale) alongside other model parameters. Two additional
  MCMC parameters (`si_shape`, `si_scale`) are sampled. Priors: shape ~
  Uniform(0.1, 10), scale ~ Uniform(0.1, 20). The existing `serial_density()`
  C++ function computes the discretized PMF from the Weibull parameters at
  each MCMC iteration.

* **Missing onset time imputation.** Infected household contacts with missing
  (`NA`) onset times are now automatically imputed during MCMC via Bayesian
  data augmentation, drawing uniformly from the follow-up window and
  accepting/rejecting via the full likelihood. An informational message
  reports the number and percentage of missing onsets. Index case onset
  times must still be non-missing.

* **SI defaults to bundled influenza serial interval.** The `SI` argument in
  `household_dynamics()`, `run_MCMC()`, and `simulate_data()` now defaults to
  `NULL`, which loads the bundled flu SI from Tsang et al. (2014). Users no
  longer need to call `data(SI)` and pass it explicitly.

## Bug fixes

* `simulate_data()` now forces single-threaded execution to avoid
  thread-unsafe `R::runif` calls in `parallelFor`. This fixes correlated
  random draws that attenuated covariate effects in simulated data.

## Other

* Added testthat test suite (47 tests across 4 files).
* `run_MCMC()` now warns when `with_rm = 1` (experimental random-effects
  branch with known per-household recording limitation).

# hhdynamics 1.1.0

## New features

* **Missing covariate imputation via Bayesian data augmentation.** Factor
  covariates with missing values (NA) are now automatically imputed during
  MCMC, jointly with onset time imputation. The MCMC draws from the full
  conditional distribution at each iteration using a uniform categorical
  proposal over factor levels. No user action is required — simply pass data
  with NAs in factor covariates. An informational message reports the number
  and percentage of missing values per covariate.

## Limitations

* Only categorical (factor) covariates are imputed. Continuous covariates
  with missing values will still produce an error.
* Interaction terms (`~sex*age`) are not supported with missing covariate data.
* The same variable cannot appear in both `inf_factor` and `sus_factor` if it
  has missing values.

# hhdynamics 1.0.1

## Breaking changes

* `household_dynamics()` now returns an S3 object of class `hhdynamics_fit`
  instead of a plain data frame. Use `summary(fit)` to get the old-style
  parameter estimates table.

* Formula arguments (`inf_factor`, `sus_factor`) now use R formula objects
  (`~sex`, `~age`) instead of character strings (`'~sex'`, `'~age'`). Use
  `NULL` instead of `'~'` when no covariates are needed.

## New features

* S3 methods for `hhdynamics_fit` objects:
    - `print()`: brief model summary (data size, MCMC settings, runtime)
    - `summary()`: parameter estimates with credible intervals
    - `coef()`: named vector of posterior means

* Full MCMC output preserved in the fit object. Access via:
    - `fit$samples`: posterior parameter samples
    - `fit$log_likelihood`: log-likelihood trace (for convergence diagnostics)
    - `fit$acceptance`: per-parameter acceptance rates
    - `fit$update_accept`: per-iteration infection-time update acceptance
    - `fit$imputed_data`: final imputed data matrix

* Input validation in `household_dynamics()` catches common errors
  (missing columns, wrong formula types, bad MCMC settings) with clear
  messages before the C++ MCMC runs.

## Bug fixes

* Fixed buffer overflow in C++ MCMC when household size exceeds 20.
  The `updateacceptrate` matrix is now sized to `max_member` instead of
  a hardcoded 20.

* `para_summary()` no longer overwrites the median with the mean, and
  no longer has plotting side effects (`layout()`, `par()`).

* `run_MCMC()` now uses `message()` instead of `print()` for runtime output.

# hhdynamics 1.0

* Initial release.
