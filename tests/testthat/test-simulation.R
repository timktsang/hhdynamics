test_that("simulate_data returns correct dimensions", {
  data(inputdata)
  data(SI)

  para <- c(1, 0.01, 0.1, 0, 0.1, 0.1, 0.1)
  sim <- simulate_data(inputdata, rep_num = 2, ~sex, ~age,
                       SI = SI, para = para, with_rm = 0)

  expect_true(is.matrix(sim))
  expect_equal(nrow(sim), 2 * 386)  # 386 households x 2 reps
})

test_that("simulate_data works without covariates", {
  data(inputdata)
  data(SI)

  para <- c(1, 0.01, 0.1, 0)
  sim <- simulate_data(inputdata, rep_num = 1, NULL, NULL,
                       SI = SI, para = para, with_rm = 0)

  expect_true(is.matrix(sim))
  expect_equal(nrow(sim), 386)
})

test_that("simulate_data runs single-threaded (RNG safety)", {
  data(inputdata)
  data(SI)

  para <- c(1, 0.01, 0.1, 0, -0.3, 0.2, -0.5)

  # Run twice with same seed — should produce identical results
  set.seed(123)
  sim1 <- simulate_data(inputdata, rep_num = 1, ~sex, ~age,
                        SI = SI, para = para, with_rm = 0)
  set.seed(123)
  sim2 <- simulate_data(inputdata, rep_num = 1, ~sex, ~age,
                        SI = SI, para = para, with_rm = 0)

  expect_identical(sim1, sim2)
})

test_that("simulate_data restores thread count", {
  prev <- RcppParallel::defaultNumThreads()
  data(inputdata)
  data(SI)

  para <- c(1, 0.01, 0.1, 0)
  simulate_data(inputdata, 1, NULL, NULL, SI = SI, para = para, with_rm = 0)

  expect_equal(RcppParallel::defaultNumThreads(), prev)
})
