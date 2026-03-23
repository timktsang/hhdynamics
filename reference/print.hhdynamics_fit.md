# Print method for hhdynamics_fit

Print method for hhdynamics_fit

## Usage

``` r
# S3 method for class 'hhdynamics_fit'
print(x, ...)
```

## Arguments

- x:

  An object of class `hhdynamics_fit`.

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
print(fit)
#> Household transmission model fit
#>   Data: 386 households, 1533 individuals
#>   MCMC: 500 post-burnin samples (1000 iterations, burnin: 500, thin: 1)
#>   Runtime: 2 seconds
#>   Parameters: community, household, sex1.0, age1.0, age2.0
#> 
#> Use summary() for estimates, coef() for posterior means.
# }
```
