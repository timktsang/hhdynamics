#' Example of input data
#'
#' This is an example of the input data used in the \code{hhdynamics} function. This data frame illustrates the format of the input data. This is a simulated data.
#' @docType data
#' @usage data(inputdata)
#' @format A example data with 8 variables, where each row represents an individual. For user's data, more predictors could be added by adding columns. For categorical variables, it should be specificed as factor, by using as.factor() function. The date is the dataset is an integer, by selecting a reference date as day 1. In the example dataset, 2008-01-01 is day 1.
#' \describe{
#'   \item{hhID}{The household id}	
#'   \item{member}{The member id in a households, use 0 to index, 1 and so on as household contact}	
#'   \item{size}{The number of individuals in the household}
#'   \item{end}{The end date of follow-up of that individual}
#'   \item{inf}{The infection status of the member. By defintion, index must be infected}	
#'   \item{onset}{The onset time of infected individual}	
#'   \item{age_group}{The age group of individual. 0: 0-19, 1: 20-64, 2: 65 or above}
#'   \item{sex}{The sex of the individual. 0: Female, 1: Male}
#' }
#' @family inputdata
"inputdata"

#' Example of serial interval distribution of influenza
#'
#' This is an example of the serial interval distribution used in the \code{hhdynamics} function. This vector specifies the format of the serial interval distribution. This is estimated from Tsang et al. Association between antibody titers and protection against influenza virus infection within households. J Infect Dis. 2014 Sep 1;210(5):684-92.
#' @docType data
#' @usage data(SI)
#' @format A vector of length 14, where element X of the vector represents a probability that the onset day is the X days after the onset day of the infector.
#' \describe{
#'   \item{}{This is the serial interval distribution. The sum of the vector elements should be 1.}
#'}
#' @family example_data
"SI"

#' Example of parameter vector for the main model
#'
#' This is an example of the parameter vector for the main model used in the \code{hhdynamics} function. This vector specifies the format of the parameter vector for the main model.
#' @docType data
#' @usage data(para)
#' @format A vector with 7 elements, where each of them is a model parameter:
#' \describe{
#'   \item{element 1}{the random effect of individual infectivity (not available in this version. Models assumed no indiviudal heterogeneity after accounting for factors affecting infectivity)}
#'   \item{element 2}{the probability of infection from the community}
#'   \item{element 3}{the probability of person-to-person transmission in households}
#'   \item{element 4}{the parameter of the relationship between the number of household contacts and transmission (default 0 in this version)}
#'   \item{element 5}{the parameter for infectivity of male (Reference group: male)}
#'   \item{element 6}{the relative susceptibility of age group 1 (Reference group: age group 0)}
#'   \item{element 7}{the relative susceptibility of age group 2 (Reference group: age group 0)}
#'}
#' @family example_data
"para"
