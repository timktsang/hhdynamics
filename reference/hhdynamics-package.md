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
#>  Probability of person-to-person transmission in households          0.058
#>                                                      sex1.0          0.013
#>                                                      age1.0         -0.153
#>                                                      age2.0         -0.376
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>        0.002       0.007                  NA               NA               NA
#>        0.030       0.082                  NA               NA               NA
#>       -0.610       0.676               1.013            0.543            1.966
#>       -0.617       0.284               0.858            0.539            1.329
#>       -0.842       0.092               0.686            0.431            1.097
# }
```
