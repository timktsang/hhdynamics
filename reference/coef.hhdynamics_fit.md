# Extract model coefficients from hhdynamics_fit

Returns named vector of posterior means for all estimated parameters.

## Usage

``` r
# S3 method for class 'hhdynamics_fit'
coef(object, ...)
```

## Arguments

- object:

  An object of class `hhdynamics_fit`.

- ...:

  Additional arguments (unused).

## Value

A named numeric vector of posterior means.

## Examples

``` r
# \donttest{
data(inputdata)
fit <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
coef(fit)
#>        re_sd    community    household   size_param       sex1.0       age1.0 
#>  1.000000000  0.004445209  0.050958903  0.000000000  0.041475617 -0.009232523 
#>       age2.0 
#> -0.278304519 
# }
```
