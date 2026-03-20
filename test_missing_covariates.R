# Test script for hhdynamics v1.1.0 — missing covariate imputation
# Run from the repo root: Rscript test_missing_covariates.R

library(hhdynamics)
data(inputdata)
data(SI)

set.seed(123)
N_ITER <- 10000
BURNIN <- 2000
THIN <- 5

cat("============================================================\n")
cat("hhdynamics v1.1.0 — Missing Covariate Imputation Test Suite\n")
cat("============================================================\n\n")

# ── Test 1: Baseline (complete data) ──────────────────────────────
cat("── Test 1: Complete data baseline ──\n")
fit_complete <- household_dynamics(inputdata, ~sex, ~age, SI,
                                   n_iteration = N_ITER, burnin = BURNIN, thinning = THIN)
s_complete <- summary(fit_complete)
cat("Complete data estimates:\n")
print(s_complete)
cat("\n")

# ── Test 2: Missing sus covariate (age, ~2% missing) ─────────────
cat("── Test 2: Missing susceptibility covariate (age) ──\n")
inputdata_age_na <- inputdata
na_idx_age <- sample(which(inputdata_age_na$member != 0), 30)
inputdata_age_na$age[na_idx_age] <- NA
cat(sprintf("  Introduced %d NAs in 'age' (%.1f%%)\n",
            sum(is.na(inputdata_age_na$age)),
            100 * sum(is.na(inputdata_age_na$age)) / nrow(inputdata_age_na)))

fit_age_na <- household_dynamics(inputdata_age_na, ~sex, ~age, SI,
                                 n_iteration = N_ITER, burnin = BURNIN, thinning = THIN)
s_age_na <- summary(fit_age_na)
cat("Missing age estimates:\n")
print(s_age_na)
cat("\n")

# ── Test 3: Missing inf covariate (sex, ~2% missing) ─────────────
cat("── Test 3: Missing infectivity covariate (sex) ──\n")
inputdata_sex_na <- inputdata
na_idx_sex <- sample(which(inputdata_sex_na$member != 0), 30)
inputdata_sex_na$sex[na_idx_sex] <- NA
cat(sprintf("  Introduced %d NAs in 'sex' (%.1f%%)\n",
            sum(is.na(inputdata_sex_na$sex)),
            100 * sum(is.na(inputdata_sex_na$sex)) / nrow(inputdata_sex_na)))

fit_sex_na <- household_dynamics(inputdata_sex_na, ~sex, ~age, SI,
                                 n_iteration = N_ITER, burnin = BURNIN, thinning = THIN)
s_sex_na <- summary(fit_sex_na)
cat("Missing sex estimates:\n")
print(s_sex_na)
cat("\n")

# ── Test 4: Both covariates missing (~5% each) ───────────────────
cat("── Test 4: Both covariates missing ──\n")
inputdata_both <- inputdata
na_idx_both_age <- sample(which(inputdata_both$member != 0), 50)
na_idx_both_sex <- sample(which(inputdata_both$member != 0), 50)
inputdata_both$age[na_idx_both_age] <- NA
inputdata_both$sex[na_idx_both_sex] <- NA
cat(sprintf("  %d NAs in 'age', %d NAs in 'sex'\n",
            sum(is.na(inputdata_both$age)), sum(is.na(inputdata_both$sex))))

fit_both <- household_dynamics(inputdata_both, ~sex, ~age, SI,
                               n_iteration = N_ITER, burnin = BURNIN, thinning = THIN)
s_both <- summary(fit_both)
cat("Both missing estimates:\n")
print(s_both)
cat("\n")

# ── Test 5: Sentinel and factor group check ───────────────────────
cat("── Test 5: Wide data sentinel check ──\n")
res <- create_wide_data(inputdata_both, ~sex, ~age)
cat(sprintf("  factor_group: [%s]\n", paste(res[[4]], collapse = ", ")))
cat(sprintf("  n_levels_vec: [%s]\n", paste(res[[5]], collapse = ", ")))
cat(sprintf("  -99 sentinels: %d\n", sum(as.matrix(res[[1]]) == -99)))
cat(sprintf("  -1 padding:    %d\n", sum(as.matrix(res[[1]]) == -1)))
cat("\n")

# ── Test 6: Continuous covariate with NA → error ──────────────────
cat("── Test 6: Continuous covariate NA → expected error ──\n")
inputdata_cont <- inputdata
inputdata_cont$cont <- rnorm(nrow(inputdata_cont))
inputdata_cont$cont[5] <- NA
tryCatch(
  household_dynamics(inputdata_cont, ~sex, ~cont, SI, n_iteration = 100, burnin = 50),
  error = function(e) cat(sprintf("  OK: %s\n", e$message))
)
cat("\n")

# ── Test 7: Interaction + missing → error ─────────────────────────
cat("── Test 7: Interaction with missing → expected error ──\n")
tryCatch(
  household_dynamics(inputdata_age_na, ~sex * age, SI = SI, n_iteration = 100, burnin = 50),
  error = function(e) cat(sprintf("  OK: %s\n", e$message))
)
cat("\n")

# ── Test 8: Shared variable with missing in both formulas → error ─
cat("── Test 8: Same variable in both formulas with NA → expected error ──\n")
inputdata_shared <- inputdata
inputdata_shared$sex[5] <- NA
tryCatch(
  household_dynamics(inputdata_shared, ~sex, ~sex, SI, n_iteration = 100, burnin = 50),
  error = function(e) cat(sprintf("  OK: %s\n", e$message))
)
cat("\n")

# ── Test 9: No covariates (backward compat) ──────────────────────
cat("── Test 9: No covariates (backward compatibility) ──\n")
fit_nocov <- household_dynamics(inputdata, SI = SI,
                                n_iteration = N_ITER, burnin = BURNIN, thinning = THIN)
cat("No-covariate estimates:\n")
print(summary(fit_nocov))
cat("\n")

# ── Comparison table ──────────────────────────────────────────────
cat("============================================================\n")
cat("Comparison: community and household probability estimates\n")
cat("============================================================\n")
get_probs <- function(s) {
  comm <- s[s$Variable == "Daily probability of infection from community", ]
  hh   <- s[s$Variable == "Probability of person-to-person transmission in households", ]
  c(comm = comm[["Point estimate"]], hh = hh[["Point estimate"]])
}
comp <- data.frame(
  scenario = c("Complete", "Age NA (2%)", "Sex NA (2%)", "Both NA (5%)", "No covariates"),
  community = c(get_probs(s_complete)["comm"],
                get_probs(s_age_na)["comm"],
                get_probs(s_sex_na)["comm"],
                get_probs(s_both)["comm"],
                get_probs(summary(fit_nocov))["comm"]),
  household = c(get_probs(s_complete)["hh"],
                get_probs(s_age_na)["hh"],
                get_probs(s_sex_na)["hh"],
                get_probs(s_both)["hh"],
                get_probs(summary(fit_nocov))["hh"])
)
rownames(comp) <- NULL
print(comp, digits = 4)
cat("\nExpected: estimates should be similar across scenarios (missing data\n")
cat("imputation should not substantially change point estimates).\n")
cat("\nDone.\n")
