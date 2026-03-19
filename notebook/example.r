## example code to fit household transmission model to the example data.

rm(list = ls())

library(devtools)
install_github("timktsang/hhdynamics")
library(hhdynamics)

data("inputdata")
data("SI")

### fit household transmission model
fit <- household_dynamics(inputdata, ~sex, ~age, SI,
  n_iteration = 15000, burnin = 5000, thinning = 1)

# Brief overview
print(fit)

# Parameter estimates table
summary(fit)

# Posterior means
coef(fit)

# Access full MCMC samples
dim(fit$samples)

# Log-likelihood trace
dim(fit$log_likelihood)
