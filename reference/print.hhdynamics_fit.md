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
print(fit)
#> Household transmission model fit
#>   Data: 386 households, 1533 individuals
#>   MCMC: 10000 post-burnin samples (15000 iterations, burnin: 5000, thin: 1)
#>   Runtime: 41 seconds
#>   Parameters: community, household, sex1.0, age1.0, age2.0
#> 
#> Use summary() for estimates, coef() for posterior means.
# }
```
