test_that("NA onset for infected contacts is accepted and imputed", {
  data(inputdata, package = "hhdynamics")
  infected_contacts <- which(inputdata$member > 0 & inputdata$inf == 1)
  inp <- inputdata
  inp$onset[infected_contacts[1:2]] <- NA
  expect_message(
    fit <- household_dynamics(inp, ~sex, ~age,
      n_iteration = 500, burnin = 100, thinning = 1),
    "missing onset"
  )
  expect_s3_class(fit, "hhdynamics_fit")
})

test_that("NA onset for uninfected contacts is fine", {
  data(inputdata, package = "hhdynamics")
  uninfected <- which(inputdata$member > 0 & inputdata$inf == 0)
  inp <- inputdata
  inp$onset[uninfected[1:3]] <- NA
  expect_no_error(
    fit <- household_dynamics(inp, n_iteration = 500, burnin = 100)
  )
})

test_that("index case onset cannot be NA", {
  data(inputdata, package = "hhdynamics")
  inp <- inputdata
  idx_cases <- which(inp$member == 0)
  inp$onset[idx_cases[1]] <- NA
  expect_error(household_dynamics(inp), "Index cases.*non-missing onset")
})

test_that("results comparable with and without missing onset", {
  skip_on_cran()
  data(inputdata, package = "hhdynamics")
  fit_full <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 5000, burnin = 2000, thinning = 1)

  infected_contacts <- which(inputdata$member > 0 & inputdata$inf == 1)
  inp <- inputdata
  inp$onset[infected_contacts[1:2]] <- NA
  fit_miss <- household_dynamics(inp, ~sex, ~age,
    n_iteration = 5000, burnin = 2000, thinning = 1)

  # Estimates should be in same ballpark (loose check)
  comm_full <- mean(fit_full$samples[, "community"])
  comm_miss <- mean(fit_miss$samples[, "community"])
  expect_true(comm_miss / comm_full > 0.3 && comm_miss / comm_full < 3.0)
})
