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
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
plot(fit)                          # diagnostics (default)

plot(fit, type = "transmission")   # transmission probability curve

# }
```
