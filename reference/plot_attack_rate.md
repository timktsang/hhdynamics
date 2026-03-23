# Plot secondary attack rates

Forest plot showing observed secondary attack rates (SAR) with Wilson
score 95% confidence intervals. Supports overall SAR, stratification by
one or more covariates, and combinations thereof. Variable names are
used as bold section headers; strata are indented below. The layout
mirrors
[`plot_covariates()`](https://timktsang.github.io/hhdynamics/reference/plot_covariates.md):
labels on the left, point estimates and CI bars in the middle, and `n/N`
counts plus SAR percentages on the right.

## Usage

``` r
plot_attack_rate(
  fit,
  by = NULL,
  include_overall = FALSE,
  labels = NULL,
  file = NULL,
  width = 8,
  height = NULL,
  xlim = NULL,
  cex = 0.85,
  ...
)
```

## Arguments

- fit:

  An object of class `hhdynamics_fit`.

- by:

  Formula, character string, or a *list* of formulas naming the
  stratification variable(s). Examples: `~age`, `"age"`,
  `list(~sex, ~age)`. Default: `NULL` (overall SAR only).

- include_overall:

  Logical. When `TRUE`, an "Overall" row is prepended even when `by` is
  specified. Default: `FALSE`.

- labels:

  Optional named list of custom display labels for variables and their
  levels. Each element is a list with `name` (display name for the
  section header) and/or `levels` (character vector of level labels in
  the same order as `sort(unique(variable))`). Names must match variable
  names in the data. Example:
  `list(age = list(name = "Age Group", levels = c("0-5", "6-17", "18+")))`.

- file:

  Optional file path for PDF output. Height is auto-calculated from the
  number of rows. Default: `NULL` (current device).

- width:

  PDF width in inches. Default: 8.

- height:

  PDF height in inches. Default: `0.45 * n_rows + 1.8`.

- xlim:

  Numeric vector of length 2 for the x-axis range (probability scale).
  Default: auto-determined from the data.

- cex:

  Character expansion factor. Default: 0.85.

- ...:

  Additional graphical parameters passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisible data frame of the estimate rows (Stratum, N_contacts,
N_infected, SAR, Lower, Upper).

## Examples

``` r
# \donttest{
data(inputdata)
fit <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 1000, burnin = 500, thinning = 1)
#> The running time is 2 seconds

# Overall only
plot_attack_rate(fit)


# Stratified by age with section header
plot_attack_rate(fit, by = ~age)


# Combined: overall + age + sex in one figure
plot_attack_rate(fit, by = list(~sex, ~age), include_overall = TRUE,
  labels = list(sex = list(name = "Sex", levels = c("Male", "Female")),
                age = list(name = "Age Group", levels = c("0-5", "6-17", "18+"))))

# }
```
