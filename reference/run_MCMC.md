# Run the MCMC for the household transmission model

Low-level function that runs the C++ MCMC sampler. Most users should use
[`household_dynamics`](https://timktsang.github.io/hhdynamics/reference/household_dynamics.md)
instead, which handles data preparation, validation, and returns an S3
object with named parameters.

## Usage

``` r
run_MCMC(
  data_w,
  SI = NULL,
  n_iteration = 15000,
  burnin = 5000,
  thinning = 1,
  n_inf,
  n_sus,
  with_rm,
  factor_group = integer(0),
  n_levels_vec = integer(0),
  estimate_SI = FALSE
)
```

## Arguments

- data_w:

  The input data, in wide format (each row is a household), as produced
  by
  [`create_wide_data`](https://timktsang.github.io/hhdynamics/reference/create_wide_data.md).

- SI:

  The mass function of the serial interval distribution. Defaults to the
  bundled influenza serial interval from Tsang et al. (2014).

- n_iteration:

  The number of iterations of the MCMC.

- burnin:

  The number of burn-in iterations to discard.

- thinning:

  The thinning interval for posterior samples.

- n_inf:

  The number of parameters affecting infectivity in the model.

- n_sus:

  The number of parameters affecting susceptibility in the model.

- with_rm:

  Indicator if the model has a random effect on individual
  infectivity (1) or not (0). **Experimental:** when `with_rm = 1`, the
  random-effects output records one value per household (index case
  only), not per individual. A warning is issued at runtime.

- factor_group:

  Integer vector mapping each dummy covariate column to its original
  factor group (from
  [`create_wide_data`](https://timktsang.github.io/hhdynamics/reference/create_wide_data.md)).

- n_levels_vec:

  Integer vector of the number of levels for each dummy column's factor
  (from
  [`create_wide_data`](https://timktsang.github.io/hhdynamics/reference/create_wide_data.md)).

- estimate_SI:

  Logical. If `TRUE`, jointly estimate Weibull shape/scale for the
  serial interval. Default is `FALSE`.

## Value

A list with 6 elements from the C++ MCMC:

1.  Posterior samples matrix (post-burnin, thinned)

2.  Log-likelihood matrix (full chain, 3 columns: total, component 1,
    component 2)

3.  Random effect samples (post-burnin). When `with_rm = 0`, returns a
    zero-variance placeholder matrix. When `with_rm = 1`
    (**experimental**), returns one random-effect value per household
    (index case only, not per individual).

4.  Acceptance rates (per-parameter, numeric vector)

5.  Infection-time update acceptance rates (iterations x household
    members)

6.  Final imputed data matrix

## Details

The MCMC uses a Metropolis-Hastings algorithm with adaptive proposal
variances. After 500 iterations, proposal standard deviations are set to
the empirical posterior standard deviation, with multiplicative tuning
based on acceptance rate (target: 20–30%). Infection times for household
contacts are jointly updated at each iteration via a data augmentation
step.

The parameter vector has the following structure:

1.  Standard deviation of random effect on infectivity (fixed at initial
    value if `with_rm = 0`)

2.  Rate of infection from community (log scale)

3.  Rate of person-to-person transmission in households (log scale)

4.  Household size parameter (currently fixed at 0)

5.  Infectivity covariate effects (`n_inf` parameters)

6.  Susceptibility covariate effects (`n_sus` parameters)

## See also

[`household_dynamics`](https://timktsang.github.io/hhdynamics/reference/household_dynamics.md)
for the high-level interface,
[`create_wide_data`](https://timktsang.github.io/hhdynamics/reference/create_wide_data.md)
for data preparation.

## Examples

``` r
# \donttest{
result_list <- create_wide_data(inputdata, ~sex, ~age)
data_w <- result_list[[1]]
n_inf <- result_list[[2]]
n_sus <- result_list[[3]]
mcmc_result <- run_MCMC(data_w,
  n_iteration = 1000, burnin = 500,
  thinning = 1, n_inf = n_inf, n_sus = n_sus, with_rm = 0)
#> The running time is 2 seconds
# }
```
