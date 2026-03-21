test_that("household_dynamics returns hhdynamics_fit with covariates", {
  data(inputdata)
  data(SI)

  fit <- household_dynamics(inputdata, ~sex, ~age, SI,
                            n_iteration = 500, burnin = 100, thinning = 1)

  expect_s3_class(fit, "hhdynamics_fit")
  expect_equal(ncol(fit$samples), 7)  # re_sd, community, household, size, sex1, age1, age2
  expect_equal(nrow(fit$samples), 400)  # 500 - 100
  expect_true(all(c("community", "household") %in% colnames(fit$samples)))
})

test_that("household_dynamics works without covariates", {
  data(inputdata)
  data(SI)

  fit <- household_dynamics(inputdata, SI = SI,
                            n_iteration = 500, burnin = 100, thinning = 1)

  expect_s3_class(fit, "hhdynamics_fit")
  expect_equal(ncol(fit$samples), 4)  # re_sd, community, household, size
})

test_that("summary returns correct structure", {
  data(inputdata)
  data(SI)

  fit <- household_dynamics(inputdata, ~sex, ~age, SI,
                            n_iteration = 500, burnin = 100, thinning = 1)
  s <- summary(fit)

  expect_true(is.data.frame(s))
  expect_equal(nrow(s), 5)  # community, household, sex1, age1, age2
  expect_true("Point estimate" %in% names(s))
  expect_true("Lower bound" %in% names(s))
  expect_true("Upper bound" %in% names(s))
})

test_that("coef returns named vector", {
  data(inputdata)
  data(SI)

  fit <- household_dynamics(inputdata, ~sex, ~age, SI,
                            n_iteration = 500, burnin = 100, thinning = 1)
  cc <- coef(fit)

  expect_true(is.numeric(cc))
  expect_true(!is.null(names(cc)))
  expect_equal(length(cc), 7)
})

test_that("print does not error", {
  data(inputdata)
  data(SI)

  fit <- household_dynamics(inputdata, SI = SI,
                            n_iteration = 500, burnin = 100, thinning = 1)
  expect_output(print(fit), "Household transmission model")
})

test_that("household_dynamics does not expose with_rm", {
  expect_false("with_rm" %in% names(formals(household_dynamics)))
})

test_that("random_effects have zero variance when with_rm=0", {
  data(inputdata)
  data(SI)

  fit <- household_dynamics(inputdata, ~sex, ~age, SI,
                            n_iteration = 200, burnin = 50, thinning = 5)

  re_var <- apply(fit$random_effects, 2, var)
  expect_true(all(re_var == 0))
})
