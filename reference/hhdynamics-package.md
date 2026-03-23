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
fit <- household_dynamics(inputdata, inf_factor = ~sex, sus_factor = ~age,
  n_iteration = 1000, burnin = 500)
#> The running time is 2 seconds
summary(fit)
#>                                                    Variable Point estimate
#>               Daily probability of infection from community          0.004
#>  Probability of person-to-person transmission in households          0.059
#>                                                      sex1.0         -0.123
#>                                                      age1.0         -0.074
#>                                                      age2.0         -0.366
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>        0.002       0.007                  NA               NA               NA
#>        0.036       0.082                  NA               NA               NA
#>       -0.733       0.363               0.884            0.480            1.438
#>       -0.597       0.371               0.929            0.550            1.449
#>       -0.866       0.146               0.693            0.421            1.157
# }
```
