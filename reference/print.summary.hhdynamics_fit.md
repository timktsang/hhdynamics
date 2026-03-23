# Print method for summary.hhdynamics_fit

Print method for summary.hhdynamics_fit

## Usage

``` r
# S3 method for class 'summary.hhdynamics_fit'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `summary.hhdynamics_fit`.

- digits:

  Number of significant digits for printing. Default is 3.

- ...:

  Additional arguments (unused).

## Value

Invisibly returns `x`.

## Examples

``` r
# \donttest{
data(inputdata)
fit <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
s <- summary(fit)
print(s, digits = 4)
#>                                                    Variable Point estimate
#>               Daily probability of infection from community         0.0046
#>  Probability of person-to-person transmission in households         0.0572
#>                                                      sex1.0        -0.1211
#>                                                      age1.0        -0.1090
#>                                                      age2.0        -0.3642
#>  Lower bound Upper bound exp(Point estimate) exp(Lower bound) exp(Upper bound)
#>       0.0021      0.0080                  NA               NA               NA
#>       0.0359      0.0824                  NA               NA               NA
#>      -0.6974      0.5291              0.8860           0.4979           1.6975
#>      -0.5111      0.3359              0.8967           0.5998           1.3991
#>      -0.8350      0.0561              0.6947           0.4339           1.0577
# }
```
