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
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
plot_transmission(fit)

# }
```
