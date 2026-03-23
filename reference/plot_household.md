# Plot household infection timeline

Visualizes the infection timeline for a single household. The index case
is shown as a filled triangle, infected contacts as filled circles at
their (imputed) onset times, and uninfected contacts as open circles
spanning their follow-up period.

## Usage

``` r
plot_household(fit, hh_id, col = NULL, ...)
```

## Arguments

- fit:

  An object of class `hhdynamics_fit`.

- hh_id:

  Household identifier to visualize.

- col:

  Colors for infected and uninfected members. Default:
  `c("firebrick", "grey60")`.

- ...:

  Additional graphical parameters.

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
plot_household(fit, hh_id = 1)

# }
```
