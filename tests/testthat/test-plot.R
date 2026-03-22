test_that("plot_diagnostics runs without error", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  pdf(NULL)
  expect_no_error(plot_diagnostics(fit))
  dev.off()
})

test_that("plot_diagnostics with specific params", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  pdf(NULL)
  expect_no_error(plot_diagnostics(fit, params = "community"))
  dev.off()
})

test_that("plot_diagnostics errors on unknown param", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  expect_error(plot_diagnostics(fit, params = "nonexistent"), "Unknown parameter")
})

test_that("plot(fit) dispatches to diagnostics", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  pdf(NULL)
  expect_no_error(plot(fit))
  dev.off()
})

test_that("plot_transmission returns data frame", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  pdf(NULL)
  result <- plot_transmission(fit)
  dev.off()
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 14)
  expect_true(all(c("Day", "Median", "Lower", "Upper") %in% names(result)))
})

test_that("plot_attack_rate returns data frame", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  pdf(NULL)
  result <- plot_attack_rate(fit, by = ~age)
  dev.off()
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 1)
})

test_that("plot_household runs without error", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata, ~sex, ~age,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  hh_ids <- unique(inputdata$hhID)
  pdf(NULL)
  expect_no_error(plot_household(fit, hh_id = hh_ids[1]))
  dev.off()
})

test_that("plot_household errors on invalid hh_id", {
  data(inputdata, package = "hhdynamics")
  fit <- household_dynamics(inputdata,
    n_iteration = 3000, burnin = 1000, thinning = 1)
  expect_error(plot_household(fit, hh_id = -999), "not found")
})
