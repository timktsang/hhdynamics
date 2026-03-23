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
#> The running time is 45 seconds
coef(fit)
#>        re_sd    community    household   size_param       sex1.0       age1.0 
#>  1.000000000  0.004292916  0.057520681  0.000000000 -0.066102155 -0.061204425 
#>       age2.0 
#> -0.307227386 
# }
```
