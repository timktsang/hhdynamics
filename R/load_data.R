#' Example of input data
#'
#' This is an example of the input data used in the \code{hhdynamics} function. This data frame illustrates the format of the input data. This is a simulated data.
#' @docType data
#' @usage data(inputdata)
#' @format A data frame with 8 variables, where each row represents an individual:
#' \describe{
#'   \item{hhID}{The household id}	
#'   \item{member}{The member id in a households, use 0 to index, 1 and so on as household contact}	
#'   \item{size}{The number of individuals in the household}
#'   \item{end}{The end date of follow-up of that individual}
#'   \item{inf}{The infection status of the member. By defintion, index must be infected}	
#'   \item{onset}{The onset time of infected individual.}	
#'   \item{age_group}{The age group of individual. 0: 0-9, 1: 10-19, 2: 20-49, 3: 50-64, 4: 65 or above }
#'   \item{sex}{The sex of the individual. 0: Female, 1: Male}
#' }
#' @family inputdata
"inputdata"

#' Example of serial interval distribution of influenza
#'
#' This is an example of the serial interval distribution used in the \code{hhdynamics} function. This vector specifies the format of the serial interval distribution. This is estimated from Tsang et al. Association between antibody titers and protection against influenza virus infection within households. J Infect Dis. 2014 Sep 1;210(5):684-92.
#' @docType data
#' @usage data(flu_activity)
#' @format A vector of length 14, where element X of the vector represents a probability that the onset day is the X days after infection.
#' \describe{
#'   \item{}{This is the serial interval distribution. The sum of the vector elements should be 1.}
#'}
#' @family example_data
"SI"

