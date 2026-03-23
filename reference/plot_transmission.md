# Plot transmission probability over time since onset

Shows the daily probability of person-to-person transmission as a
function of days since the infector's symptom onset. The serial interval
distribution shapes this curve. The median and 95% credible interval are
computed from the posterior samples.

## Usage

``` r
plot_transmission(fit, hh_size = NULL, col = "steelblue", ...)
```

## Arguments

- fit:

  An object of class `hhdynamics_fit`.

- hh_size:

  Reference household size. Default: median from data.

- col:

  Color for the median line and credible interval polygon. Default:
  `"steelblue"`.

- ...:

  Additional graphical parameters passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisible data frame with columns Day, Median, Lower, Upper.

## Details

When the model was fitted with `estimate_SI = TRUE`, the uncertainty
band incorporates serial interval uncertainty (via the Weibull
shape/scale posterior).

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
plot_transmission(fit)

# }
```
