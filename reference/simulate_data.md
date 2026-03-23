# Simulate household transmission data

Generates synthetic datasets from the household transmission model for
validation, power analysis, or posterior predictive checks.

## Usage

``` r
simulate_data(
  input,
  rep_num,
  inf_factor = NULL,
  sus_factor = NULL,
  SI = NULL,
  para,
  with_rm
)
```

## Arguments

- input:

  The dataset in long format (same structure as for
  [`household_dynamics`](https://timktsang.github.io/hhdynamics/reference/household_dynamics.md)).

- rep_num:

  The number of replications of the input dataset, to increase the
  sample size.

- inf_factor:

  Formula for factors affecting infectivity (e.g. `~sex`). Use `NULL`
  for no factors. Default is `NULL`.

- sus_factor:

  Formula for factors affecting susceptibility (e.g. `~age`). Use `NULL`
  for no factors. Default is `NULL`.

- SI:

  The mass function of the serial interval distribution. Defaults to the
  bundled influenza serial interval from Tsang et al. (2014).

- para:

  The parameter vector, matching the structure from
  [`coef.hhdynamics_fit`](https://timktsang.github.io/hhdynamics/reference/coef.hhdynamics_fit.md): (1)
  random effect SD, (2) community rate, (3) household rate, (4) size
  parameter, (5+) covariate effects.

- with_rm:

  Indicator if the model has a random effect on individual
  infectivity (1) or not (0).

## Value

A simulated dataset in wide format (one row per household) based on the
input parameter vectors.

## Details

The simulation uses the same household structure (sizes, follow-up
periods, covariate values) as the input data. The `rep_num` parameter
replicates the household structure to increase sample size. Infection
outcomes and onset times are simulated from the model given the
parameter vector.

The output is in wide format (one row per household), matching the
internal representation used by the C++ backend. Use this with
[`run_MCMC()`](https://timktsang.github.io/hhdynamics/reference/run_MCMC.md)
or
[`household_dynamics()`](https://timktsang.github.io/hhdynamics/reference/household_dynamics.md)
to verify model recovery.

## See also

[`household_dynamics`](https://timktsang.github.io/hhdynamics/reference/household_dynamics.md)
for fitting the model,
[`coef.hhdynamics_fit`](https://timktsang.github.io/hhdynamics/reference/coef.hhdynamics_fit.md)
for extracting parameter estimates to use as simulation inputs.

## Examples

``` r
# \donttest{
data(inputdata)
data(SI)
para <- c(1, 0.01, 0.1, 0, 0.1, 0.1, 0.1)
simulated <- simulate_data(inputdata, 10, ~sex, ~age,
  SI = SI, para = para, with_rm = 0)
# }
```
