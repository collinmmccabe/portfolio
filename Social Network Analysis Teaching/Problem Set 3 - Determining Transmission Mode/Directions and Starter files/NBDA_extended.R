#------------------------------------------------------------#
#                      NBDA_extended.r                       #
# - Compares models of mixed social-asocial & exclusively    #
# - asocial transmission through a predefined network.       #
#                                                            #
#             [Adapted from Franz & Nunn 2009]               #
#                                             ~Collin McCabe #
#------------------------------------------------------------#


#----------------------#
#      Functions       #
#----------------------#

# function that calculates the log-likelihood for the non-social learning model for a specific parameterization; this function is used as input for an optimization routine that is provided in R [M Franz]
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

# function that calculates the log-likelihood for a model of non-social and social transmission, for a specific parameterization; this function is used as input for an optimization routine that is provided in R [M Franz]
log_lik_sam <- function(par, social_network, diffusion_data, time_max)
{
  alr <- par[1]
  tau <- par[2]
  
  if (alr >= 0 && alr <= 1 && tau >= 0 )
  { 
    log_lik <- 0
    group_size <- dim(social_network)[1]
    
    for (time in 1:time_max)
    {  
      skilled <- diffusion_data < time    # individuals that were already infected in the previous time step
      learning_rates <- 1 - exp(-1 * tau * colSums( ((skilled %*% matrix(1,1,group_size)) * social_network) ) )  # these are the probabilities to become infected socially for each individual
      
      for (i in 1:group_size)
      {
        if (diffusion_data[i] == time)  #if this individual learned in the current time step
          log_lik <- log_lik + log(1 - ((1 - alr) * (1 - learning_rates[i]))) # this gives the log-likelihood that the current individual was infected by social or non-social transmission                                                              
        if (diffusion_data[i] > time)  # if this individual was not infected and did not become infected in this time step
          log_lik <- log_lik + log((1 - alr) * (1 - learning_rates[i])) # this gives the log-likelihood that the current individual did not become infected socially and not by non-social contact                                                              
      }
    }
  }
  else
    log_lik <- -Inf
  
  return(log_lik)
}

# function that calculates AIC values and akaike weights for the social and the non-social learning model, given the observed diffusion data and a network [M Franz]
extended_nbda <- function(social_network, diffusion_data, tau_i, alr_i, time_max)
{
  
  # maximum likelihood estimation for the non-social transmission model
  alr_ml <- optimize(f = log_lik_am, lower = 0, upper = 1, maximum = TRUE, tol = 0.001, diffusion_data  = diffusion_data, time_max = time_max)
  max_log_lik_a <- alr_ml$objective
  
  # maximum likelihood estimation for the model of social and non-social transmission
  tau_alr_ml <- optim(par = c(alr_i, tau_i), fn = log_lik_sam, control = list(fnscale = - 1), social_network = social_network, diffusion_data  = diffusion_data, time_max = time_max)
  max_log_lik_sa <- tau_alr_ml$value
  
  # calculate AIC values and Akaike weights
  results = data.frame(AIC = c(0,0), Akaike_weight = c(0,0), row.names =c("social and non-social transmission model", "pure non-social transmission model") )
  
  results$AIC[1] =  -2 * max_log_lik_sa + 2 * 2
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

tau_i <- 0        # initial value of tau that is used for estimation of tau in the model of social and non-social transmission
alr_i <- 0.5      # initial value of non-social transmission rate that is used for estimation of non-social transmission rate in the model of social and non-social transmission
time_max <- max(diff_data, na.rm=TRUE)  # this parameterization has to be changed if the trait did not spread through the whole group

diff_data[is.na(diff_data)] <- time_max + 1 # control for NAs in dataset [C McCabe]

################
# Run analysis #
################

res <- extended_nbda(social_network = SN, diffusion_data = diff_data, time_max = time_max, tau_i = tau_i, alr_i = alr_i)

# display results [C McCabe]
cat("\n"); print(res); cat("\n")