test_that("estimate_SI = FALSE (default) has no SI params", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 500, burnin = 100, thinning = 1)
  expect_false("si_shape" %in% colnames(fit$samples))
  expect_false("si_scale" %in% colnames(fit$samples))
})

test_that("estimate_SI = TRUE adds si_shape and si_scale", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 500, burnin = 100, thinning = 1,
    estimate_SI = TRUE)
  expect_true("si_shape" %in% colnames(fit$samples))
  expect_true("si_scale" %in% colnames(fit$samples))
  # SI params should be positive
  expect_true(all(fit$samples[, "si_shape"] > 0))
  expect_true(all(fit$samples[, "si_scale"] > 0))
})

test_that("estimate_SI = TRUE works without covariates", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata,
    n_iteration = 500, burnin = 100, thinning = 1,
    estimate_SI = TRUE)
  expect_true("si_shape" %in% colnames(fit$samples))
  expect_equal(ncol(fit$samples), 6)  # 4 base + 2 SI
})

test_that("estimate_SI = TRUE skips SI validation", {
  data(inputdata, package = "hhdynamics")
  # Bad SI should not error when estimate_SI = TRUE
  expect_no_error(
    household_dynamics(inputdata, SI = rep(2, 14),
      n_iteration = 500, burnin = 100, thinning = 1,
      estimate_SI = TRUE)
  )
})

test_that("summary works with SI params (no exp columns)", {
  skip_on_cran()
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 500, burnin = 100, thinning = 1,
    estimate_SI = TRUE)
  s <- summary(fit)
  # SI rows should have NA in exp columns
  si_rows <- which(s$Variable %in% c("si_shape", "si_scale"))
  expect_true(all(is.na(s[si_rows, "exp(Point estimate)"])))
})

test_that("SI defaults to bundled data when NULL", {
  data(inputdata, package = "hhdynamics")
  # Should work without specifying SI
  fit <- household_dynamics(inputdata, n_iteration = 500, burnin = 100)
  expect_s3_class(fit, "hhdynamics_fit")
})
