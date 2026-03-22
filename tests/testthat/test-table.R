test_that("table_parameters returns expected columns", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  tab <- table_parameters(fit)
  expect_s3_class(tab, "data.frame")
  expect_true(all(c("Parameter", "Mean", "Median", "Lower", "Upper", "ESS", "Acceptance") %in% names(tab)))
  # Should have community, household, 3 covariates = 5 rows (no re_sd, no size_param)
  expect_equal(nrow(tab), 5)
})

test_that("table_parameters works without covariates", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  tab <- table_parameters(fit)
  expect_equal(nrow(tab), 2)  # community + household only
})

test_that("table_covariates returns correct structure", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  tab <- table_covariates(fit)
  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 3)  # sex1, age1, age2
  expect_true(all(c("Covariate", "Type", "Estimate", "exp_Estimate") %in% names(tab)))
  expect_equal(tab$Type[1], "Infectivity")
  expect_equal(tab$Type[2], "Susceptibility")
})

test_that("table_covariates returns empty df when no covariates", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  tab <- table_covariates(fit)
  expect_equal(nrow(tab), 0)
})

test_that("table_attack_rates overall", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  tab <- table_attack_rates(fit)
  expect_equal(nrow(tab), 1)
  expect_equal(tab$Stratum, "Overall")
  expect_true(tab$SAR > 0 && tab$SAR < 1)
  expect_true(tab$Lower <= tab$SAR && tab$SAR <= tab$Upper)
})

test_that("table_attack_rates by covariate", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  tab <- table_attack_rates(fit, by = ~age)
  expect_equal(nrow(tab), 3)  # 3 age levels
  expect_true(all(tab$SAR >= 0 & tab$SAR <= 1))
})

test_that("table_attack_rates errors on missing variable", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  expect_error(table_attack_rates(fit, by = ~nonexistent), "not found")
})

test_that("table_attack_rates handles NA covariates correctly", {
  data(inputdata, package = "hhdynamics")
  d <- inputdata
  d$age[c(5, 10, 15)] <- NA
  fit <- household_dynamics(d, ~sex, ~age,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  tab <- table_attack_rates(fit, by = ~age)
  # Should only have non-NA strata (no NA row)
  expect_false(any(is.na(tab$Stratum)))
  # No NA in infected counts
  expect_false(any(is.na(tab$N_infected)))
  # Contact total should not include NA rows
  expect_true(all(tab$SAR >= 0 & tab$SAR <= 1))
})

test_that("ESS is near n for iid draws", {
  set.seed(42)
  x <- rnorm(2000)
  ess <- hhdynamics:::.effective_sample_size(x)
  # For iid draws, ESS should be close to n (at least 80%)
  expect_gt(ess, 0.8 * 2000)
})
