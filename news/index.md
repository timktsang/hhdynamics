# Changelog

## hhdynamics 1.3.1

### Bug fixes

- Replaced C-style variable length arrays (VLAs) with
  `std::vector<double>` in C++ MCMC code to comply with ISO C++
  standard. Fixes CRAN check WARNING on Windows (GCC 14) and
  installation failure on Debian Linux.

- Removed obsolete `CXX_STD = CXX11` from Makevars/Makevars.win and
  `SystemRequirements: C++11` from DESCRIPTION. C++11 support has been
  removed in R-devel.

- Added LAPACK/BLAS linking flags to `src/Makevars` (previously only in
  `Makevars.win`), fixing `undefined symbol: dpotrf_` installation error
  on Debian Linux.

## hhdynamics 1.3.0

### New features

- **Forest plot for covariate effects.** `plot_covariates(fit)` produces
  a publication-ready forest plot grouped by Susceptibility and
  Infectiousness sections, with reference level labels, alternating row
  shading, and a dashed null line at RR = 1. Optional `file` argument
  saves a PDF with auto-sized dimensions based on the number of
  covariate rows. Custom variable and level labels can be supplied via
  the `labels` argument.

- **Summary tables.** Three new exported functions return clean data
  frames: `table_parameters(fit)` (posterior mean/median/CrI with ESS
  and acceptance), `table_covariates(fit)` (covariate effects on the log
  and exponentiated scale), and `table_attack_rates(fit, by)` (secondary
  attack rates with Wilson CIs, optionally stratified by a covariate).

### Bug fixes

- Fixed
  [`table_attack_rates()`](https://timktsang.github.io/hhdynamics/reference/table_attack_rates.md)
  NA handling: rows with `NA` in the stratifying variable are now
  correctly excluded from strata and counts (previously `x == NA`
  returned `NA` rather than `FALSE`, poisoning subset sizes).

- Fixed `.effective_sample_size()` lag indexing: Geyer (1992) pairing
  now correctly starts at lag 1 (index 2 of R’s
  [`acf()`](https://rdrr.io/r/stats/acf.html) output) rather than lag 0,
  which was inflating the integrated autocorrelation time by ~2 and
  causing ESS to be systematically underestimated by roughly half.

## hhdynamics 1.2.0

### New features

- **Optional serial interval estimation.** Set `estimate_SI = TRUE` in
  [`household_dynamics()`](https://timktsang.github.io/hhdynamics/reference/household_dynamics.md)
  to jointly estimate the serial interval distribution as a
  Weibull(shape, scale) alongside other model parameters. Two additional
  MCMC parameters (`si_shape`, `si_scale`) are sampled. Priors: shape ~
  Uniform(0.1, 10), scale ~ Uniform(0.1, 20). The existing
  `serial_density()` C++ function computes the discretized PMF from the
  Weibull parameters at each MCMC iteration.

- **Missing onset time imputation.** Infected household contacts with
  missing (`NA`) onset times are now automatically imputed during MCMC
  via Bayesian data augmentation, drawing uniformly from the follow-up
  window and accepting/rejecting via the full likelihood. An
  informational message reports the number and percentage of missing
  onsets. Index case onset times must still be non-missing.

- **SI defaults to bundled influenza serial interval.** The `SI`
  argument in
  [`household_dynamics()`](https://timktsang.github.io/hhdynamics/reference/household_dynamics.md),
  [`run_MCMC()`](https://timktsang.github.io/hhdynamics/reference/run_MCMC.md),
  and
  [`simulate_data()`](https://timktsang.github.io/hhdynamics/reference/simulate_data.md)
  now defaults to `NULL`, which loads the bundled flu SI from Tsang et
  al. (2014). Users no longer need to call `data(SI)` and pass it
  explicitly.

### Bug fixes

- [`simulate_data()`](https://timktsang.github.io/hhdynamics/reference/simulate_data.md)
  now forces single-threaded execution to avoid thread-unsafe `R::runif`
  calls in `parallelFor`. This fixes correlated random draws that
  attenuated covariate effects in simulated data.

### Other

- Added testthat test suite (47 tests across 4 files).
- [`run_MCMC()`](https://timktsang.github.io/hhdynamics/reference/run_MCMC.md)
  now warns when `with_rm = 1` (experimental random-effects branch with
  known per-household recording limitation).

## hhdynamics 1.1.0

### New features

- **Missing covariate imputation via Bayesian data augmentation.**
  Factor covariates with missing values (NA) are now automatically
  imputed during MCMC, jointly with onset time imputation. The MCMC
  draws from the full conditional distribution at each iteration using a
  uniform categorical proposal over factor levels. No user action is
  required — simply pass data with NAs in factor covariates. An
  informational message reports the number and percentage of missing
  values per covariate.

### Limitations

- Only categorical (factor) covariates are imputed. Continuous
  covariates with missing values will still produce an error.
- Interaction terms (`~sex*age`) are not supported with missing
  covariate data.
- The same variable cannot appear in both `inf_factor` and `sus_factor`
  if it has missing values.

## hhdynamics 1.0.1

### Breaking changes

- [`household_dynamics()`](https://timktsang.github.io/hhdynamics/reference/household_dynamics.md)
  now returns an S3 object of class `hhdynamics_fit` instead of a plain
  data frame. Use `summary(fit)` to get the old-style parameter
  estimates table.

- Formula arguments (`inf_factor`, `sus_factor`) now use R formula
  objects (`~sex`, `~age`) instead of character strings (`'~sex'`,
  `'~age'`). Use `NULL` instead of `'~'` when no covariates are needed.

### New features

- S3 methods for `hhdynamics_fit` objects:
  - [`print()`](https://rdrr.io/r/base/print.html): brief model summary
    (data size, MCMC settings, runtime)
  - [`summary()`](https://rdrr.io/r/base/summary.html): parameter
    estimates with credible intervals
  - [`coef()`](https://rdrr.io/r/stats/coef.html): named vector of
    posterior means
- Full MCMC output preserved in the fit object. Access via:
  - `fit$samples`: posterior parameter samples
  - `fit$log_likelihood`: log-likelihood trace (for convergence
    diagnostics)
  - `fit$acceptance`: per-parameter acceptance rates
  - `fit$update_accept`: per-iteration infection-time update acceptance
  - `fit$imputed_data`: final imputed data matrix
- Input validation in
  [`household_dynamics()`](https://timktsang.github.io/hhdynamics/reference/household_dynamics.md)
  catches common errors (missing columns, wrong formula types, bad MCMC
  settings) with clear messages before the C++ MCMC runs.

### Bug fixes

- Fixed buffer overflow in C++ MCMC when household size exceeds 20. The
  `updateacceptrate` matrix is now sized to `max_member` instead of a
  hardcoded 20.

- `para_summary()` no longer overwrites the median with the mean, and no
  longer has plotting side effects
  ([`layout()`](https://rdrr.io/r/graphics/layout.html),
  [`par()`](https://rdrr.io/r/graphics/par.html)).

- [`run_MCMC()`](https://timktsang.github.io/hhdynamics/reference/run_MCMC.md)
  now uses [`message()`](https://rdrr.io/r/base/message.html) instead of
  [`print()`](https://rdrr.io/r/base/print.html) for runtime output.

## hhdynamics 1.0

- Initial release.
