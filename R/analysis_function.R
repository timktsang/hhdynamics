#' Run the MCMC for the household transmission
#'
#' This function runs the MCMC for the household transmission model. It is used in the main function.
#' @param data_w The data for running MCMC, in dataframe format. It should be in the same format as the data in the package. It includes: 1) age_group (0: children, 1: adults, 2: older adults), 2) start_time: start of follow-up, 3) end_time: end of follow-up, 4) time1: date for first serum collection, 5) time2: date for second serum collection, 6) time3: date for third serum collection, 7) HAI_titer_1: HAI titer for first serum collection, 8) HAI_titer_2: HAI titer for second serum collection, 9) HAI_titer_3: HAI titer for third serum collection.
#' @param SI The data for influenza activity used in the inference. The row number should match with the date in the input data.
#' @param n_iteration The number of iterations of the MCMC.
#' @param burnin The iteration for burn-in for MCMC.
#' @param thinning The number of thinning in MCMC.
#' @param n_inf The number of parameters affecting infectivity in the model.
#' @param n_sus The number of parameters affecting susceptibility in the model.
#' @param with_rm Indicator if the model has a random effect on individual infectivity or not.
#' @return A matrix that stores the posterior samples for the model parameter.
#' @examples 
#' mcmc_result <- run_MCMC(data_w,SI,n_iteration = 15000,burnin = 5000,thinning = 1,n_inf,n_sus,with_rm)
#' @export
run_MCMC <- function(data_w,SI,n_iteration = 15000,burnin = 5000,thinning = 1,n_inf,n_sus,with_rm){
  
keep_iteration <- burnin + 1:((n_iteration - burnin)/thinning)*thinning
#### put to the MCMC
## first need to create parameter vector
# para 1. random effect on infectivity
# para 2. from community
# para 3. from household
# para 4. the household size parameter
# para about the input variables

para <-  c(1,0.01,0.1,0,rep(0.1,n_inf),rep(0.1,n_sus))

sigma <- c(1,rep(0.1,length(para)-1))
move <- rep(1,length(para))
move[1] <- with_rm
move[4] <- 0
sep1 <- 5
sep2 <- n_inf+n_sus+3 # 3 for onset time and the inf status, and the random effect

aaaaa1 <- Sys.time()
tt <- mcmc(as.matrix(data_w),SI,n_iteration,burnin,thinning,para,move,sigma,n_inf,n_sus,with_rm,sep1,sep2)
aaaaa2 <- Sys.time()
print(aaaaa2-aaaaa1)

return(tt[[1]][keep_iteration,])
}

#################################################
## function to compute the mcmc output
## internal use, not output
para_summary <- function(mcmc,a,b,print){
  y <- matrix(NA,ncol(mcmc),4)
  for (i in 1:ncol(mcmc)){
    y[i,1:3] <- quantile(mcmc[,i],c(0.5,0.025,0.975),na.rm=T)
    y[i,4] <- sum(diff(mcmc[,i])!=0)/nrow(mcmc)
    y[i,1] <- mean(mcmc[,i],na.rm=T)
  }
  layout(matrix(1:(a*b),nrow=a,byrow=T))
  par(mar=c(2,4,1,1))
  if (print==1){
    for (i in 1:ncol(mcmc)){
      plot(mcmc[,i],type="l")
    }
  }
  return(y)
}

#################################################
#' Create a wide format of household data
#'
#' This function creates a wide format of the household data (each row is a household), based on the long format of data (each row is an individual), to fit the household transmission model.
#' @param input The input data, in long format.
#' @param inf_factor Factors affecting infectivity. Use the format of '~ factor1 + factor2'.
#' @param sus_factor Factors affecting susceptibility. Use the format of '~ factor1 + factor2'.
#' @return A data frame in the wide format.
#' @examples 
#' wide_data <- create_wide_data(inputdata,'~sex','~age')
#' @export
create_wide_data <- function(input,inf_factor,sus_factor){
  
  # first keep the essential
  data_a <- input[,c("hhID","member","inf","onset")]
  n_sus <- 0
  n_inf <- 0
  # get the infectivity
  data_a <- cbind(data_a,0)
  if (!(missing(inf_factor))){
    design_matrix_inf <- model.matrix(as.formula(inf_factor),input)
    n_inf <- ncol(design_matrix_inf)-1
    data_a <- cbind(data_a,design_matrix_inf[,-1,drop=F])
  }
  # get the susceptibility
  if (!missing(sus_factor)){
    design_matrix_sus <- model.matrix(as.formula(sus_factor),input)
    n_sus <- ncol(design_matrix_sus)-1
    data_a <- cbind(data_a,design_matrix_sus[,-1,drop=F])
  }
  
  
  data_w <- reshape(data_a,direction='wide',idvar='hhID',timevar = 'member')
  
  data_w <- merge(input[input$member==0,c("hhID","size","size","onset",'end')],data_w,by='hhID')
  
  data_w[is.na(data_w)] <- -1
  
  return(list(data_w,n_inf,n_sus))
}


#################################################
#' Fitting household transmission model to the data, based on Markov chain Monte Carlo (MCMC).
#'
#' The main function to fit the household transmission model to the data, to estimate the probability of infection from the community, probability of person-to-person transmission in households, and factors affecting susceptibility and infectivity.
#' @param input The input data, in long format (each row is an individual). The data must include: 1) hhID: the household id, 2) member: the member id in a household, using 0 to index, 1 and so on as household contact, 3) size: the number of individuals in the households, 4) end: the end date of follow-up for that individual, 5) inf: the infection status of the member. By definition, the index must be infected. 6) onset: the onset time of the individual. The data must also include the factors affecting infectivity or susceptibility that may be explored.
#' @param inf_factor Factors affecting infectivity. Use the format of '~ factor1 + factor2'.
#' @param sus_factor Factors affecting susceptibility. Use the format of '~ factor1 + factor2'.
#' @param SI The mass function of the serial interval distribution. It should have a length equal to 14.
#' @param n_iteration The number of iterations of the MCMC.
#' @param burnin The iteration for burn-in for MCMC.
#' @param thinning The number of thinning in MCMC.
#' @return A data frame that stores the estimates from a fitted MCMC.
#' @examples 
#' fitted_result <- household_dynamics(inputdata,'~sex','~age',SI,15000,5000,1,0)
#' @export
household_dynamics <- function(input,inf_factor,sus_factor,SI,n_iteration = 15000,burnin = 5000,thinning = 1){
  # later may add the model with random effect on infectivity.
  with_rm <- 0 
  # first keep the essential
  result_list <- create_wide_data(input,inf_factor,sus_factor)
 
  data_w <- result_list[[1]]
  n_inf <- result_list[[2]]
  n_sus <- result_list[[3]]
  
testing2 <- run_MCMC(data_w,SI,n_iteration,burnin,thinning,n_inf,n_sus,with_rm)

z1 <- para_summary(testing2,4,3,0)

output <- data.frame( matrix(NA,nrow(z1),7) )
names(output) <- c("Variable","Point estimate","Lower bound","Upper bound","exp(Point estimate)","exp(Lower bound)","exp(Upper bound)")
z1[2:3,] <- 1-exp(-z1[2:3,])
z1 <- cbind(z1,exp(z1))
output[] <- z1[,c(4,1:3,5:7)]
output[,1] <- c("Standard deviation of random effect",
                "Daily probability of infection from community",
                "Probability of person-to-person transmission in households",
                "Parameter of the relationship between transmission and number of household members",names(data_w)[5+3+(1:(n_sus+n_inf))])
output[1:4,5:7] <- NA
output <- output[-4,]
if (with_rm == 0){
output <- output[-1,]  
}

output_print <- output
output_print[,-1] <- round(  output_print[,-1],3 )
print(output_print)

return(output)
}

#################################################
#' Simulation of the dataset based on the input dataset and a household transmission model.
#'
#' This function simulates the dataset for validation or other purposes. 
#' @param input The dataset in long format.
#' @param rep_num The number of replications of the input dataset, to increase the sample size. (needs modification later)
#' @param inf_factor Factors affecting infectivity. Use the format of '~ factor1 + factor2'.
#' @param sus_factor Factors affecting susceptibility. Use the format of '~ factor1 + factor2'.
#' @param para The parameter vector for the model parameters, in the following format: 1) the random effect of individual infectivity, 2) the probability of infection from the community, 3) the probability of person-to-person transmission in households, 4) the parameter of the relationship between the number of household contacts and transmission, 5 or more: the parameters of relative infectivity or susceptibility.
#' @param with_rm Indicator if the model has a random effect on individual infectivity or not.
#' @return A simulated data based on the input parameter vectors, with the format equal to the input data.
#' @examples 
#' a1 <- simulate_data(inputdata,10,'~sex','~age',SI,c(1,0.005,0.05,0,0.1,0.1,0.1,0.1,0.1),0)
#' @export
simulate_data <- function(input,rep_num,inf_factor,sus_factor,SI,para,with_rm){
 
  create_sim_data <- create_wide_data(input,inf_factor,sus_factor)
  simdata <- create_sim_data[[1]]
  
  n_inf <- create_sim_data[[2]]
  n_sus <- create_sim_data[[3]]

  c1 <- sim_data(as.matrix(simdata[rep(1:nrow(simdata),rep_num),]),SI,para,n_inf,n_sus,with_rm,5,n_inf+n_sus+3)
  
  return(c1[[1]])
}
