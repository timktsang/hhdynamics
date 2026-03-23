# Plot method for hhdynamics_fit objects

Plot method for hhdynamics_fit objects

## Usage

``` r
# S3 method for class 'hhdynamics_fit'
plot(x, type = "diagnostics", ...)
```

## Arguments

- x:

  An object of class `hhdynamics_fit`.

- type:

  Type of plot: `"diagnostics"` (default), `"transmission"`,
  `"attack_rate"`, or `"covariates"`.

- ...:

  Additional arguments passed to the underlying plot function.

## Value

Invisible NULL.

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
plot(fit)                          # diagnostics (default)

plot(fit, type = "transmission")   # transmission probability curve

# }
```
