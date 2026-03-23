# Fit a household transmission model via MCMC

The main function to fit the household transmission model to data.
Estimates the daily probability of infection from the community, the
probability of person-to-person transmission within households, and
effects of covariates on infectivity and susceptibility.

## Usage

``` r
household_dynamics(
  input,
  inf_factor = NULL,
  sus_factor = NULL,
  SI = NULL,
  n_iteration = 15000,
  burnin = 5000,
  thinning = 1,
  estimate_SI = FALSE
)
```

## Arguments

- input:

  The input data, in long format (each row is an individual). Required
  columns:

  hhID

  :   Household identifier.

  member

  :   Member index (0 = index case, 1+ = contacts).

  size

  :   Number of individuals in the household.

  end

  :   End date of follow-up for that individual.

  inf

  :   Infection status (1 = infected, 0 = not). Index cases must have
      `inf = 1`.

  onset

  :   Onset time of symptoms.

- inf_factor:

  Formula for factors affecting infectivity (e.g. `~sex` or
  `~sex + age`). Use `NULL` for no factors. Default is `NULL`.

- sus_factor:

  Formula for factors affecting susceptibility (e.g. `~age`). Use `NULL`
  for no factors. Default is `NULL`.

- SI:

  The mass function of the serial interval distribution. Must be a
  numeric vector of length 14 summing to approximately 1. Defaults to
  the bundled influenza serial interval from Tsang et al. (2014). Not
  used when `estimate_SI = TRUE` (SI is estimated from data via Weibull
  parameterization).

- n_iteration:

  Total number of MCMC iterations. Default is 15000.

- burnin:

  Number of initial iterations to discard. Default is 5000.

- thinning:

  Thinning interval for posterior samples. Default is 1.

- estimate_SI:

  Logical. If `TRUE`, jointly estimate the serial interval distribution
  as a Weibull(shape, scale) alongside other model parameters. Two
  additional parameters (`si_shape`, `si_scale`) are added to the MCMC.
  Priors: shape ~ Uniform(0.1, 10), scale ~ Uniform(0.1, 20). Default is
  `FALSE`.

## Value

An object of class
[`print.hhdynamics_fit`](https://timktsang.github.io/hhdynamics/reference/print.hhdynamics_fit.md)`{hhdynamics_fit}`.
Use [`summary()`](https://rdrr.io/r/base/summary.html) to get parameter
estimates, [`print()`](https://rdrr.io/r/base/print.html) for a brief
overview, and [`coef()`](https://rdrr.io/r/stats/coef.html) for
posterior means. When `estimate_SI = TRUE`, the output includes
`si_shape` and `si_scale` parameters.

## Details

The model assumes that each household contact can be infected either
from the community (at a constant daily rate) or from an infected
household member (with probability governed by the serial interval
distribution). Tertiary transmission within households is accounted for.
Infection times for non-index cases are treated as latent variables and
imputed via data augmentation during MCMC.

Covariate effects on infectivity and susceptibility enter
multiplicatively on the log scale. The
[`summary()`](https://rdrr.io/r/base/summary.html) method reports
exponentiated estimates for interpretation as relative risks.

**Missing covariate imputation:** Factor covariates with missing values
(`NA`) are automatically imputed during MCMC via Bayesian data
augmentation, using a uniform categorical prior over factor levels. Only
factor covariates are supported; continuous covariates with `NA` will
produce an error. Interaction terms with missing data are not supported.

The returned `hhdynamics_fit` object stores the full MCMC output,
enabling custom convergence diagnostics and post-processing. Key fields:

- `$samples`:

  Posterior parameter samples (post-burnin, thinned). Columns named by
  parameter.

- `$log_likelihood`:

  Log-likelihood trace for convergence assessment (full chain).

- `$acceptance`:

  Per-parameter acceptance rates from the Metropolis-Hastings sampler.

- `$imputed_data`:

  Final imputed dataset (wide format) with augmented infection times.

## See also

[`simulate_data`](https://timktsang.github.io/hhdynamics/reference/simulate_data.md)
for simulating from the model,
[`create_wide_data`](https://timktsang.github.io/hhdynamics/reference/create_wide_data.md)
for data preparation,
[`run_MCMC`](https://timktsang.github.io/hhdynamics/reference/run_MCMC.md)
for the low-level MCMC interface.

## Examples

``` r
# \donttest{
data(inputdata)

# Fit with default flu SI
fit <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
print(fit)
#> Household transmission model fit
#>   Data: 386 households, 1533 individuals
#>   MCMC: 500 post-burnin samples (1000 iterations, burnin: 500, thin: 1)
#>   Runtime: 2 seconds
#>   Parameters: community, household, sex1.0, age1.0, age2.0
#> 
#> Use summary() for estimates, coef() for posterior means.
summary(fit)
#>                                                    Variable Point estimate
#>               Daily probability of infection from community          0.004
#>  Probability of person-to-person transmission in households          0.057
#>                                                      sex1.0         -0.087
#>                                                      age1.0         -0.053
#>                                                      age2.0         -0.328
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>        0.003       0.007                  NA               NA               NA
#>        0.036       0.082                  NA               NA               NA
#>       -0.656       0.442               0.917            0.519            1.555
#>       -0.491       0.412               0.949            0.612            1.510
#>       -0.753       0.121               0.721            0.471            1.129
coef(fit)
#>       re_sd   community   household  size_param      sex1.0      age1.0 
#>  1.00000000  0.00422478  0.05847233  0.00000000 -0.08667671 -0.05260881 
#>      age2.0 
#> -0.32780242 
# }
```
