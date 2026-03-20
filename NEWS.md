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
