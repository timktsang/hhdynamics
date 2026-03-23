# Limit threads to 2 for CRAN compliance
RcppParallel::setThreadOptions(numThreads = 2L)

# Shared fixtures: fit once, reuse across plot/table/diagnostics tests
data(inputdata, package = "hhdynamics")

.fit_cov <- household_dynamics(inputdata, ~sex, ~age,
  n_iteration = 500, burnin = 100, thinning = 1)

.fit_nocov <- household_dynamics(inputdata,
  n_iteration = 500, burnin = 100, thinning = 1)
