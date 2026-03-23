test_that("plot_diagnostics runs without error", {
  pdf(NULL)
  expect_no_error(plot_diagnostics(.fit_cov))
  dev.off()
})

test_that("plot_diagnostics with specific params", {
  pdf(NULL)
  expect_no_error(plot_diagnostics(.fit_cov, params = "community"))
  dev.off()
})

test_that("plot_diagnostics errors on unknown param", {
  expect_error(plot_diagnostics(.fit_nocov, params = "nonexistent"), "Unknown parameter")
})

test_that("plot(fit) dispatches to diagnostics", {
  pdf(NULL)
  expect_no_error(plot(.fit_nocov))
  dev.off()
})

test_that("plot_transmission returns data frame", {
  pdf(NULL)
  result <- plot_transmission(.fit_cov)
  dev.off()
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 14)
  expect_true(all(c("Day", "Median", "Lower", "Upper") %in% names(result)))
})

test_that("plot_attack_rate returns data frame", {
  pdf(NULL)
  result <- plot_attack_rate(.fit_cov, by = ~age)
  dev.off()
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 1)
})

test_that("plot_household runs without error", {
  data(inputdata, package = "hhdynamics")
  hh_ids <- unique(inputdata$hhID)
  pdf(NULL)
  expect_no_error(plot_household(.fit_cov, hh_id = hh_ids[1]))
  dev.off()
})

test_that("plot_household errors on invalid hh_id", {
  expect_error(plot_household(.fit_nocov, hh_id = -999), "not found")
})
