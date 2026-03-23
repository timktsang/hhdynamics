# Forest plot of covariate effects

Produces a forest plot showing estimated relative risks for covariate
effects on susceptibility and infectiousness. Covariates are grouped by
variable with bold headers, reference categories labeled, alternating
row shading, and estimate text with credible intervals on the right.

## Usage

``` r
plot_covariates(
  fit,
  probs = c(0.025, 0.975),
  labels = NULL,
  file = NULL,
  width = 11,
  height = NULL,
  xlim = NULL,
  xlab_left = "Lower Risk",
  xlab_right = "Higher Risk",
  cex = 0.85,
  ...
)
```

## Arguments

- fit:

  An object of class `hhdynamics_fit`.

- probs:

  Numeric vector of length 2 for credible interval bounds. Default:
  `c(0.025, 0.975)` (95% CrI).

- labels:

  Optional named list of custom labels for covariates. Each element is a
  list with `name` (display name for the variable header) and `levels`
  (character vector of level labels, including the reference level
  first). Names must match variable names in the formula. Example:
  `list(sex = list(name = "Sex", levels = c("Male", "Female")))`.

- file:

  Optional file path for PDF output. When provided, a PDF is created
  with auto-calculated width and height based on the number of rows.
  Default: `NULL` (plot to current device).

- width:

  PDF width in inches. Default: 11. Only used when `file` is not `NULL`.

- height:

  PDF height in inches. Default: auto-calculated as
  `0.45 * n_rows + 1.8`. Only used when `file` is not `NULL`.

- xlim:

  Numeric vector of length 2 for the x-axis range on the natural scale.
  Default: auto-determined from the credible intervals.

- xlab_left:

  Label for the left direction arrow. Default: `"Lower Risk"`.

- xlab_right:

  Label for the right direction arrow. Default: `"Higher Risk"`.

- cex:

  Character expansion factor. Default: 0.85.

- ...:

  Additional graphical parameters passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisible NULL. Called for its side effect of producing a plot.

## Details

When `file` is provided, the plot is saved to a PDF with dimensions
automatically calculated from the number of covariate rows. When `file`
is `NULL`, the plot is drawn to the current graphics device.

## Examples

``` r
# \donttest{
data(inputdata)
fit <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds
plot_covariates(fit)


# Save to PDF with auto-sized dimensions
plot_covariates(fit, file = tempfile(fileext = ".pdf"),
  labels = list(sex = list(name = "Sex", levels = c("Male", "Female")),
                age = list(name = "Age Group", levels = c("0-5", "6-17", "18+"))))
# }
```
