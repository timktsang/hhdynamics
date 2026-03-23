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
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
plot_household(fit, hh_id = 1)

# }
```
