test_that("missing factor covariate in sus_factor is accepted and imputed", {
  data(inputdata)
  data(SI)

  inputdata_na <- inputdata
  na_idx <- sample(which(inputdata_na$member != 0), 30)
  inputdata_na$age[na_idx] <- NA

  expect_message(
    fit <- household_dynamics(inputdata_na, ~sex, ~age, SI,
                              n_iteration = 500, burnin = 100, thinning = 1),
    "will be imputed"
  )
  expect_s3_class(fit, "hhdynamics_fit")
})

test_that("missing factor covariate in inf_factor is accepted and imputed", {
  data(inputdata)
  data(SI)

  inputdata_na <- inputdata
  na_idx <- sample(which(inputdata_na$member != 0), 30)
  inputdata_na$sex[na_idx] <- NA

  expect_message(
    fit <- household_dynamics(inputdata_na, ~sex, ~age, SI,
                              n_iteration = 500, burnin = 100, thinning = 1),
    "will be imputed"
  )
  expect_s3_class(fit, "hhdynamics_fit")
})

test_that("missing continuous covariate raises error", {
  data(inputdata)
  data(SI)

  inputdata_cont <- inputdata
  inputdata_cont$cont <- rnorm(nrow(inputdata_cont))
  inputdata_cont$cont[5] <- NA

  expect_error(
    household_dynamics(inputdata_cont, ~sex, ~cont, SI, n_iteration = 100, burnin = 50),
    "Only factor covariates"
  )
})

test_that("interaction with missing data raises error", {
  data(inputdata)
  data(SI)

  inputdata_na <- inputdata
  na_idx <- sample(which(inputdata_na$member != 0), 10)
  inputdata_na$age[na_idx] <- NA

  expect_error(
    household_dynamics(inputdata_na, ~sex * age, SI = SI, n_iteration = 100, burnin = 50),
    "Interaction terms"
  )
})

test_that("shared variable with missing in both formulas raises error", {
  data(inputdata)
  data(SI)

  inputdata_shared <- inputdata
  inputdata_shared$sex[5] <- NA

  expect_error(
    household_dynamics(inputdata_shared, ~sex, ~sex, SI, n_iteration = 100, burnin = 50),
    "appears in both"
  )
})

test_that("wide data sentinels are correct with missing covariates", {
  data(inputdata)

  inputdata_na <- inputdata
  inputdata_na$age[sample(which(inputdata_na$member != 0), 20)] <- NA
  inputdata_na$sex[sample(which(inputdata_na$member != 0), 20)] <- NA

  res <- create_wide_data(inputdata_na, ~sex, ~age)
  wide <- as.matrix(res[[1]])

  # Should have -99 sentinels for missing covariates
  expect_true(sum(wide == -99) > 0)
  # Should have -1 padding for absent household members
  expect_true(sum(wide == -1) > 0)
  # factor_group and n_levels_vec should be populated
  expect_equal(length(res[[4]]), 3)  # sex(1) + age(2)
  expect_equal(length(res[[5]]), 3)
})
