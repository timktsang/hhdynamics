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
#>  Probability of person-to-person transmission in households          0.055
#>                                                      sex1.0         -0.088
#>                                                      age1.0         -0.037
#>                                                      age2.0         -0.268
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>        0.002       0.007                  NA               NA               NA
#>        0.035       0.078                  NA               NA               NA
#>       -0.725       0.465               0.916            0.484            1.592
#>       -0.531       0.378               0.964            0.588            1.459
#>       -0.894       0.170               0.765            0.409            1.185
# }
```
