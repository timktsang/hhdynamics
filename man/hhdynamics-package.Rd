\name{hhdynamics-package}
\alias{hhdynamics-package}
\alias{hhdynamics}
\docType{package}
\title{
   Household Transmission Modeling for household transmission (Case-Ascertain) Studies
}
\description{
  The `hhdynamics` package provides a model for analyzing household transmission dynamics, specifically focusing on fitting models to case-ascertained household transmission studies. It allows users to assess the risk of infections among household contacts based on various parameters. The package also provides functions for simulating household transmission dynamics and for fitting models to data from case-ascertained household transmission studies.
}
\details{
This package includes functions to simulate household transmission data, fit models using Markov Chain Monte Carlo (MCMC) methods, and generate summary statistics for the analysis. 
Main function is:
   - `household_dynamics()`: high-level function to get estimates from modeling of household transmission studies. 
Key functions include:
  - `simulate_data()`: Generates synthetic data for household transmission studies.
  - `run_MCMC()`: Fits the transmission model using MCMC.
  - `create_wide_data()`: Prepares data for model fitting.
  - `SI()`: Calculates secondary infection probabilities.'
}
\author{
Tim Tsang, timkltsang@gmail.com.

Maintainer: Tim Tsang <timkltsang@gmail.com>
}
\references{
To cite package **hhdynamics** in publications use:

Tsang TK, Cauchemez S, Perera RA, Freeman G, Fang VJ, Ip DK, Leung GM,
Malik Peiris JS, Cowling BJ. (2014). Association between antibody titers
and protection against influenza virus infection within households. J
Infect Dis. 2014 Sep 1;210(5):684-92
}
\keyword{ package }
\seealso{
  README.md
}
\examples{
  \dontrun{
## Load in the example data and the serial interval distribution.
data("inputdata")
data("SI")

###### run the MCMC to estimate parameter of the model
###### in actual analysis, number of iteration is 15000, burnin is 5000, and thinning is 1
fit_hh_model <- household_dynamics(inputdata, inf_factor = '~sex', sus_factor = '~age',SI,n_iteration =  15000,burnin = 5000, thinning =  1)
  }
}
