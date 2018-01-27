#------------------------------------------------------------#
#                      NBDA_basic.r                          #
# - Compares models of both exclusively social & exclusively #
# - asocial transmission through a predefined network.       #
#                                                            #
#             [Adapted from Franz & Nunn 2009]               #
#                                             ~Collin McCabe #
#------------------------------------------------------------#


#----------------------#
#      Functions       #
#----------------------#

# function that calculates the log-likelihood for the social transmission model for a specific parameterization; this function is used as input for an optimization routine that is provided in R [M Franz]
log_lik_sm <- function(tau, social_network, diffusion_data, time_max)
{
  log_lik <- 0
  group_size <- dim(social_network)[1]
  for (time in 1:time_max)
  {
    skilled <- diffusion_data < time    # individuals that were already infected in the previous time step
    learning_rates <- 1 - exp(-1 * tau * colSums( ((skilled %*% matrix(1,1,group_size)) * social_network) ) )
    
    log_lik <- log_lik + sum( log(learning_rates * (diffusion_data == time) + (diffusion_data != time)) )    # individuals that were infected in this time step, all transmission rates of individuals that did not become infected in this time step are automatically set to 1, thus they do not affect the log-likelihood
    log_lik <- log_lik + sum( log(( 1 - learning_rates) * (diffusion_data > time) + (diffusion_data <= time)))    # individuals that did not yet become infected, all transmission rates of individuals that are already infected are automatically set to 1, thus they do not affect the log-likelihood
  }
  
  return(log_lik)
}

# function that calculates the log-likelihood for the non-social transmission model for a specific parameterization this function is used as input for an optimization routine that is provided in R [M Franz]
log_lik_am <- function(alr, diffusion_data, time_max)
{
  log_lik <- 0
  for (time in 1:time_max)
  {      
    log_lik <- log_lik + sum(diffusion_data == time) * log(alr)     # individuals that were infected in this time step
    log_lik <- log_lik + sum(diffusion_data > time) * log(1 - alr)    # individuals that did not yet become infected
  }
  
  return(log_lik)
}

# function that calculates AIC values and Akaike weights given the observed diffusion data and a network [M Franz]
nbda <- function(social_network, diffusion_data, tau_max, time_max)
{
  
  # maximum likelihood estimation for the non-social transmission model
  alr_ml <- optimize(f = log_lik_am, lower = 0, upper = 1, maximum = TRUE, tol = 0.001, diffusion_data  = diffusion_data, time_max = time_max)
  max_log_lik_a <- alr_ml$objective
  
  # maximum likelihood estimation for the social transmission model
  tau_ml <- optimize(f = log_lik_sm, lower = 0, upper = tau_max, maximum = TRUE, tol = 0.001, social_network = social_network, diffusion_data  = diffusion_data, time_max = time_max)
  max_log_lik_s <- tau_ml$objective
  
  # calculate AIC values and Akaike weights
  results = data.frame(AIC = c(0,0), Akaike_weight = c(0,0), row.names =c("pure social transmission model", "pure non-social transmission model") )
  
  results$AIC[1] =  -2 * max_log_lik_s + 2
  results$AIC[2] =  -2 * max_log_lik_a + 2
  
  results$Akaike_weight[1] = exp(-0.5*(results$AIC[1] - min(results$AIC))) / ( exp(-0.5*(results$AIC[1] - min(results$AIC))) + exp(-0.5*(results$AIC[2] - min(results$AIC))) )
  results$Akaike_weight[2] = 1  - results$Akaike_weight[1]
  
  return(results)
}

# various user input prompt functions [C McCabe]
UserNetworkInput <- function() 
{
  readline("\nPlease type the filename of the sociomatrix dataset: ")
}

UserSpreadInput <- function() 
{
  readline("Please type the filename of the disease spread dataset: ")
}


#-------------------#
#   Analysis Code   #
#     [M Franz]     #
#-------------------#

###################
# Load data files #
###################

# read social network data
netfp <- UserNetworkInput()
SN <- read.csv(netfp, header=FALSE)

# read diffusion data
spreadfp <- UserSpreadInput()
diff_data <- read.csv(spreadfp, header=FALSE)

# transform this data to a column vector
diff_data <- t(diff_data)

##################
# Set parameters #
##################

tau_max <- 1  # maximum value for parameter tau in the social transmission model
time_max <- max(diff_data, na.rm=TRUE)  # this parameterization has to be changed if the trait did not spread through the whole group

diff_data[is.na(diff_data)] <- time_max + 1 # control for NAs in dataset [C McCabe]

#################
# Run analysis #
#################

options(warn=-1)
res <- nbda(social_network = SN, diffusion_data = diff_data, tau_max = tau_max, time_max = time_max)

# display results [C McCabe]
cat("\n"); print(res); cat("\n")