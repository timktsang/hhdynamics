# Table of secondary attack rates

Computes observed secondary attack rates (SAR) from the data, optionally
stratified by a covariate. Confidence intervals use the Wilson score
method.

## Usage

``` r
table_attack_rates(fit, by = NULL)
```

## Arguments

- fit:

  An object of class `hhdynamics_fit`.

- by:

  Formula or character string naming the stratification variable (e.g.
  `~age` or `"age"`). Default: no stratification.

## Value

A data frame with columns: Stratum, N_contacts, N_infected, SAR, Lower,
Upper.

## Examples

``` r
# \donttest{
data(inputdata)
fit <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 15000, burnin = 5000, thinning = 1)
#> Iteration: 1000
#> Iteration: 2000
#> Iteration: 3000
#> Iteration: 4000
#> Iteration: 5000
#> Iteration: 6000
#> Iteration: 7000
#> Iteration: 8000
#> Iteration: 9000
#> Iteration: 10000
#> Iteration: 11000
#> Iteration: 12000
#> Iteration: 13000
#> Iteration: 14000
#> The running time is 41 seconds
table_attack_rates(fit)
#>       Stratum N_contacts N_infected        SAR      Lower     Upper
#> lower Overall       1147         92 0.08020924 0.06585537 0.0973656
table_attack_rates(fit, by = ~age)
#>   Stratum N_contacts N_infected        SAR      Lower      Upper
#> 1       0        402         35 0.08706468 0.06326714 0.11867943
#> 2       1        346         30 0.08670520 0.06140856 0.12107827
#> 3       2        399         27 0.06766917 0.04691982 0.09666386
# }
```
