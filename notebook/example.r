## example code to fit household transmission model to the example data.

rm(list = ls())

library(devtools)
install_github("timktsang/hhdynamics")
library(hhdynamics)

data("inputdata")
data("SI")

### fit household transmisison model
fit_hh_model <- household_dynamics(inputdata,'~sex','~age',SI,n_iteration =  15000,burnin = 5000, thinning =  1)

