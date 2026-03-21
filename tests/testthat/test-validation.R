test_that("SI validation rejects bad inputs", {
  data(inputdata)
  data(SI)

  # Wrong sum
 expect_error(household_dynamics(inputdata, SI = rep(2, 14), n_iteration = 100, burnin = 50),
               "sum to approximately 1")

  # Negative values
  expect_error(household_dynamics(inputdata, SI = c(-0.1, rep(1/13, 13)), n_iteration = 100, burnin = 50),
               "non-negative")

  # Wrong length
  expect_error(household_dynamics(inputdata, SI = rep(0.1, 10), n_iteration = 100, burnin = 50),
               "length 14")

  # Not numeric
  expect_error(household_dynamics(inputdata, SI = rep("a", 14), n_iteration = 100, burnin = 50),
               "numeric vector")

  # All zeros
  expect_error(household_dynamics(inputdata, SI = rep(0, 14), n_iteration = 100, burnin = 50),
               "sum to approximately 1")
})

test_that("SI validation accepts valid inputs", {
  data(inputdata)
  data(SI)

  # Standard SI
  fit <- household_dynamics(inputdata, SI = SI, n_iteration = 100, burnin = 50)
  expect_s3_class(fit, "hhdynamics_fit")

  # Uniform SI
  fit2 <- household_dynamics(inputdata, SI = rep(1/14, 14), n_iteration = 100, burnin = 50)
  expect_s3_class(fit2, "hhdynamics_fit")
})

test_that("input validation catches missing columns", {
  data(SI)
  bad_data <- data.frame(x = 1:10)
  expect_error(household_dynamics(bad_data, SI = SI), "missing required columns")
})

test_that("input validation catches non-data.frame", {
  data(SI)
  expect_error(household_dynamics("not a df", SI = SI), "must be a data frame")
})

test_that("input validation catches bad formula variables", {
  data(inputdata)
  data(SI)
  expect_error(household_dynamics(inputdata, ~nonexistent, SI = SI, n_iteration = 100, burnin = 50),
               "not found in input")
})

test_that("input validation catches string formulas", {
  data(inputdata)
  data(SI)
  expect_error(household_dynamics(inputdata, "~sex", SI = SI, n_iteration = 100, burnin = 50),
               "must be a formula")
})

test_that("MCMC parameter validation works", {
  data(inputdata)
  data(SI)
  expect_error(household_dynamics(inputdata, SI = SI, n_iteration = 50, burnin = 100),
               "greater than")
  expect_error(household_dynamics(inputdata, SI = SI, n_iteration = -1, burnin = 50),
               "positive integer")
})
