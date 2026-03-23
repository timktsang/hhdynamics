# Create a wide format of household data

Reshapes long-format household data (one row per individual) into wide
format (one row per household) for the C++ MCMC backend. Also constructs
design matrices for infectivity and susceptibility covariates.

## Usage

``` r
create_wide_data(input, inf_factor, sus_factor)
```

## Arguments

- input:

  The input data, in long format. Must contain columns: hhID, member,
  inf, onset, size, end.

- inf_factor:

  Formula for factors affecting infectivity (e.g. `~sex` or
  `~sex + age`). Use `NULL` for no factors.

- sus_factor:

  Formula for factors affecting susceptibility (e.g. `~age`). Use `NULL`
  for no factors.

## Value

A list with 5 elements: (1) data frame in wide format, (2) number of
infectivity parameters, (3) number of susceptibility parameters, (4)
integer vector mapping each dummy column to its factor group, (5)
integer vector of factor levels per dummy column.

## Details

The wide-format output has the following column structure per household:

- Columns 1–5: hhID, size, size (duplicate), index onset, end date

- Columns 6+: per-member blocks of (inf, onset, random effect
  placeholder, infectivity covariates, susceptibility covariates)

Missing values (members not present in a household) are filled with -1,
which the C++ code treats as a sentinel. Missing factor covariates are
coded as -99 (a distinct sentinel) and will be imputed during MCMC.

## See also

[`household_dynamics`](https://timktsang.github.io/hhdynamics/reference/household_dynamics.md)
for the high-level interface,
[`run_MCMC`](https://timktsang.github.io/hhdynamics/reference/run_MCMC.md)
for the low-level MCMC function.

## Examples

``` r
wide_data <- create_wide_data(inputdata, ~sex, ~age)
```
