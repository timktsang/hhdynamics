# hhdynamics v1.0.1 — Refactor & Documentation Update Report

**Date:** 2026-03-19

## Summary

Refactored the hhdynamics R package (v1.0 → v1.0.1) to return a proper S3 class from the main fitting function, preserve full MCMC output, fix a C++ buffer overflow bug, modernize the formula interface, add input validation, and improve documentation with a vignette and detailed roxygen.

## Changes

### Step 1: C++ household size bug fix

**File:** `src/mcmc_function_hh.cpp:738`

The `updateacceptrate` matrix was hardcoded to 20 columns, causing a buffer overflow when any household had more than 20 members. Changed to use `max_member`, which is already computed at line 682 as `max(data1(_, 1))`.

```cpp
// Before
NumericMatrix updateacceptrate(mcmc_n, 20);

// After
NumericMatrix updateacceptrate(mcmc_n, max_member);
```

### Step 2: `run_MCMC()` returns full C++ output

**File:** `R/analysis_function.R`

Previously `run_MCMC()` discarded 5 of 6 C++ outputs, returning only the posterior samples matrix. Now returns the complete list:

| Element | Content |
|---------|---------|
| `[[1]]` | Posterior samples (post-burnin, thinned) |
| `[[2]]` | Log-likelihood trace (full chain, 3 columns) |
| `[[3]]` | Random effect samples |
| `[[4]]` | Per-parameter acceptance rates |
| `[[5]]` | Infection-time update acceptance (iterations × members) |
| `[[6]]` | Final imputed data matrix |

Also switched from `print()` to `message()` for runtime output, and attached elapsed time as an attribute.

### Step 3: Formula syntax modernized

**Files:** `R/analysis_function.R`, `notebook/example.r`

| Before | After |
|--------|-------|
| `'~sex'` | `~sex` |
| `'~age'` | `~age` |
| `'~'` (no covariates) | `NULL` |

Removed `as.formula()` wrappers in `create_wide_data()`. The old string syntax now raises a clear validation error.

### Step 4: S3 class `hhdynamics_fit`

**New file:** `R/hhdynamics_fit.R`

Created an S3 class with constructor and three methods:

| Method | Behavior |
|--------|----------|
| `print.hhdynamics_fit()` | 5-line summary: data size, MCMC settings, runtime, parameter list |
| `summary.hhdynamics_fit()` | Parameter estimates table with posterior means, 95% CIs, exp() for covariates — backward-compatible with old `household_dynamics()` output |
| `coef.hhdynamics_fit()` | Named vector of posterior means |

The `hhdynamics_fit` object stores 20 fields including all MCMC outputs, model metadata (formulas, parameter names, transforms), and run information (iterations, burnin, thinning, elapsed time, data dimensions).

### Step 5: `household_dynamics()` refactored

**File:** `R/analysis_function.R`

- Now returns `hhdynamics_fit` instead of a plain data frame
- `with_rm` exposed as an explicit parameter (default 0, previously hardcoded)
- Parameter names auto-generated from formulas + data column names
- Parameter transforms tracked (`prob` for community/household, `exp` for covariates, `none` for others)
- Old `para_summary()` call and data frame construction moved to `summary.hhdynamics_fit()`

### Step 6: Input validation

**New file:** `R/validation.R`

`validate_inputs()` checks before any C++ computation:

| Check | Error message |
|-------|---------------|
| `input` is a data.frame | `'input' must be a data frame.` |
| Required columns present | `'input' is missing required columns: ...` |
| Index cases infected | `All index cases (member == 0) must be infected (inf == 1).` |
| Formulas are formula objects or NULL | `'inf_factor' must be a formula (e.g. ~sex) or NULL.` |
| Formula variables exist in data | `Variables in inf_factor not found in input: ...` |
| SI is numeric, length 14 | `'SI' must have length 14, got N.` |
| n_iteration > burnin | `'n_iteration' must be greater than 'burnin'.` |
| with_rm is 0 or 1 | `'with_rm' must be 0 or 1.` |
| Large household size warning | `Note: maximum household size is N.` (message, not error) |

### Step 7: `para_summary()` cleaned up

**File:** `R/analysis_function.R`

- Fixed bug where median (line 58) was computed then immediately overwritten by mean (line 60) — now computes mean directly
- Removed `layout()` and `par()` side effects (trace plots will move to `plot.hhdynamics_fit()` in Phase 2)
- Returns a named data frame (`mean`, `lower`, `upper`, `acceptance`) instead of unnamed matrix

### Step 8: Documentation & packaging

| File | Change |
|------|--------|
| `DESCRIPTION` | Version 1.0 → 1.0.1, added `Suggests: knitr, rmarkdown`, `VignetteBuilder: knitr` |
| `R/hhdynamics-package.R` | Removed unused `graphics::layout`, `graphics::par`, `stats::as.formula` imports. Added `stats::sd`. |
| `NAMESPACE` | Regenerated — added S3 method registrations, removed graphics imports |
| `NEWS.md` | **New** — documents all breaking changes, features, and bug fixes |
| `vignettes/hhdynamics-intro.Rmd` | **New** — getting-started vignette covering data format, fitting, inspecting results, MCMC access, validation, simulation |
| `.Rbuildignore` | Added `Rplots.pdf` exclusion |
| All roxygen blocks | Added `@details`, `@seealso` cross-references, structured `@param` descriptions, richer `@examples` |

## Files modified

| File | Type | Steps |
|------|------|-------|
| `src/mcmc_function_hh.cpp` | Modified | 1 |
| `R/analysis_function.R` | Rewritten | 2, 3, 5, 7, 8 |
| `R/hhdynamics_fit.R` | **New** | 4 |
| `R/validation.R` | **New** | 6 |
| `R/hhdynamics-package.R` | Modified | 8 |
| `DESCRIPTION` | Modified | 8 |
| `NAMESPACE` | Regenerated | 8 |
| `NEWS.md` | **New** | 8 |
| `vignettes/hhdynamics-intro.Rmd` | **New** | 8 |
| `.Rbuildignore` | Modified | 8 |
| `notebook/example.r` | Modified | 3 |
| `man/*.Rd` | Regenerated | 8 |

## Validation

### R CMD check

```
R --no-environ CMD build .   → hhdynamics_1.0.1.tar.gz (clean)
R --no-environ CMD check --as-cran hhdynamics_1.0.1.tar.gz → 0E / 0W / 2N
```

NOTEs (both standard/acceptable):
1. "New submission" — expected for CRAN incoming feasibility
2. "unable to verify current time" — network/clock issue, harmless

### MCMC smoke test

Ran `household_dynamics(inputdata, ~sex, ~age, SI, n_iteration = 5000, burnin = 1000)` with built-in example data.

**print() output:**
```
Household transmission model fit
  Data: 386 households, 1533 individuals
  MCMC: 4000 post-burnin samples (5000 iterations, burnin: 1000, thin: 1)
  Runtime: 11 seconds
  Parameters: community, household, sex1.0, age1.0, age2.0
```

**summary() output:**
```
                                                   Variable Point estimate Lower bound Upper bound
              Daily probability of infection from community          0.004       0.002       0.007
 Probability of person-to-person transmission in households          0.057       0.034       0.084
                                                     sex1.0         -0.081      -0.733       0.488
                                                     age1.0         -0.065      -0.537       0.412
                                                     age2.0         -0.312      -0.831       0.170
```

**Dimension checks:**

| Field | Dimensions | Expected |
|-------|-----------|----------|
| `fit$samples` | 4000 × 7 | (n_iteration - burnin) × n_params |
| `fit$log_likelihood` | 5000 × 3 | n_iteration × 3 |
| `fit$acceptance` | length 7 | n_params |
| `fit$update_accept` | 5000 × 7 | n_iteration × max_member |
| `fit$imputed_data` | 386 × 47 | n_households × wide_cols |

**coef() = colMeans(fit$samples):** `TRUE`

**NULL formula path (no covariates):** Works correctly, summary shows 2 rows (community + household).

**Validation errors (all 5 caught):**
- `~nonexistent` → `Variables in inf_factor not found in input: nonexistent`
- `SI = c(0.1, 0.9)` → `'SI' must have length 14, got 2.`
- `input = "not a df"` → `'input' must be a data frame.`
- `inf_factor = "~sex"` → `'inf_factor' must be a formula (e.g. ~sex) or NULL.`
- `burnin > n_iteration` → `'n_iteration' must be greater than 'burnin'.`

## What was NOT changed

- **C++ MCMC algorithm** — only the `updateacceptrate` dimension fix; no changes to sampling logic
- **`create_wide_data()` internals** — reshape logic, NA → -1 sentinels, duplicate size column (all coupled to C++ expectations)
- **`rmrecord` column 7 hardcoding** — latent bug for random effects with certain covariate counts, documented and deferred
- **No plot/table functions** — that's Phase 2, after this foundation is solid

## Next steps

- **Phase 2:** Add `plot.hhdynamics_fit()` (trace plots, posterior densities, acceptance diagnostics)
- **Phase 2:** Add table output functions (formatted parameter tables for manuscripts)
- **Phase 1 (optional):** CRAN infra — testthat test suite, LICENSE file, CI via GitHub Actions
- Push v1.0.1 to GitHub
