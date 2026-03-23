# hhdynamics: Fitting Household Transmission Model to Estimate Household Transmission Dynamics of Influenza

A Bayesian household transmission model to estimate household
transmission dynamics, with accounting for infection from community and
tertiary cases.

## See also

Useful links:

- <https://github.com/timktsang/hhdynamics>

- Report bugs at <https://github.com/timktsang/hhdynamics/issues>

## Author

**Maintainer**: Tim Tsang <timkltsang@gmail.com>
([ORCID](https://orcid.org/0000-0001-5037-6776))

## Examples

``` r
# \donttest{
# Fit a household transmission model
data(inputdata)
fit <- household_dynamics(inputdata, inf_factor = ~sex, sus_factor = ~age)
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
#> The running time is 45 seconds
summary(fit)
#>                                                    Variable Point estimate
#>               Daily probability of infection from community          0.004
#>  Probability of person-to-person transmission in households          0.056
#>                                                      sex1.0         -0.073
#>                                                      age1.0         -0.057
#>                                                      age2.0         -0.315
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>        0.002       0.007                  NA               NA               NA
#>        0.033       0.082                  NA               NA               NA
#>       -0.713       0.514               0.930            0.490            1.672
#>       -0.532       0.425               0.944            0.588            1.530
#>       -0.829       0.174               0.730            0.436            1.190
# }
```
