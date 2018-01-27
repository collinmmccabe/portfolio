#------------------------------------------------------------#
#                      simulation.r                          #
# - Simulates the spread of a parasite through the HEB 1333  #
# - class, either socially or non-socially.                  #
#                                             ~Collin McCabe #
#------------------------------------------------------------#


#----------------------#
#      Functions       #
#----------------------#

## various user input prompt functions...
UserRnoughtInput <- function() 
{
  readline("What is the observed R0 of the parasite? ")
}

UserInfectInput <- function() 
{
  readline("What is the average observed infectious period (in hours) of the disease? ")
}

UserSpreadtypeInput <- function() 
{
  readline("What type of parasite transmission would you like to simulate?\n- input 1 for social trasnmission, or 0 for non-social transmission: ")
}

UserFilenameInput <- function() 
{
  readline("\nWhat would you like to name the output file for the spread of this parasite (must conform to <filename>.csv)? ")
}

## Install or load required packages [E Otarola-Castillo]
is_installed <- function(mypkg) { is.element(mypkg, installed.packages()[,1]) }
if(!is_installed("statnet"))
{
  install.packages("statnet")
}
library("statnet")


#-------------------#
#  Simulation Code  #
#-------------------#
 
options(warn=-1)    # suppress warnings that might arise from bogus inputs, etc.

cat("\nThis R script will run simulations of epidemics through the HEB 1333 population through either social or non-social (environmental) transmission.  You may choose any values of R0 and infectious period:\n\n")    # intro prompt

## repeat prompt with explanation for R0 if user provides bogus input
repeat {
  Rnought <- as.numeric(UserRnoughtInput())
  if(is.na(Rnought)){
    cat("\nSorry, the R0 value that you entered was not valid, try again -\n\n")
  } else {
    break
  }
}
Rnought <- round(Rnought)

## repeat prompt with explanation for infectious period if user provides bogus input
repeat {
  Inf.period <- as.numeric(UserInfectInput())
  if(is.na(Inf.period)){
    cat("\nSorry, the infectious period value that you entered was not valid, try again -\n\n")
  } else {
    break
  }
}
Inf.period <- round(Inf.period)

## repeat prompt with explanation for model type if user provides bogus input
repeat {
  model.type <- UserSpreadtypeInput()
  if(model.type != 1 && model.type != 0){
    cat("\nSorry, the model type that you entered was not valid, try again -\n\n")
  } else {
    break
  }
}

N <- 38    # population size of HEB 1333

patient.zero <- sample(1:N, 1)    # randomly assign a patient zero from population
time <- 0    # initialize a time counter

## initialize an infection time vector for later use in NBDA
infection.times <- array(NA, dim=N)
infection.times[patient.zero] <- 0

## initialize a matrix to be iteratively bound in each pass of the model (0=susceptible, 2=infected, 5=recovered)
inf.status <- array(0, dim=c(N, 1))
inf.status[patient.zero] <- 2

p.inf <- Rnought/Inf.period    # sets probablility of infection at each infectious encounter

## initialize a dataframe for keeping track of SIR modeling [Adapted from Bellam et al. 2012]
dat <- data.frame(time = time, sus = NA, inf = NA, rec = NA)
dat[1,2:4] <- c((N-1), 1, 0)

## model code for running social transmission model
if(model.type == 1) {
  associations <- read.csv("undirHEBV.csv", header=FALSE)    # import social connections from external file
  ## run model until no more infectious individuals are present in population (0mod5==0, 5mod5==0, 2mod5!=3)
  while(sum(inf.status[,(time+1)] %% 5) > 0)
  {
    tmp.inf.status <- inf.status[,(time+1)]    # make a temporary copy of the state of the population at previous time step for updating new infection statuses
    tmp.inf.status.orig <- tmp.inf.status    # make another temporary copy that will remain unchanged for reference to current state of population
    tmp.dat <- dat[(time+1),]; tmp.dat[,1] <- time    # make a temporary copy of a dataframe row for SIR model updating
    for(i in 1:N)    # for all individuals in population
    {
      ## if individual is infected and has reached the end of the infectious period...
      if(!is.na(infection.times[i]) && (infection.times[i] + Inf.period) == (time+1))
      {
        tmp.inf.status[i] <- 5    # individual becomes recovered
        tmp.dat[,3] <- tmp.dat[,3] - 1; tmp.dat[,4] <- tmp.dat[,4] + 1    # numbers of I & R updated
      ## if individual is infectious still and randomly is chosen to infect someone...
      } else if(tmp.inf.status.orig[i] == 2 && rbinom(1, 1, p.inf) == 1)
      {
        all.ties <- c(associations$V2[which(associations$V1 == i)], associations$V1[which(associations$V2 == i)])    # check to see who the individual has ties to...
        ## code to correct for limitation of sample() for a vector of length 1
        if(length(all.ties) == 1) { 
          rand.tie <- all.ties
        } else {
          rand.tie <- sample(all.ties, 1)
        }
        ## if the randomly chosen connection is susceptible, then infect
        if(tmp.inf.status[rand.tie] == 0) 
        {
          tmp.inf.status[rand.tie] <- 2    # connection becomes infected
          tmp.dat[,2] <- tmp.dat[,2] - 1; tmp.dat[,3] <- tmp.dat[,3] + 1    # numbers of S & I updated
          infection.times[rand.tie] <- time + 1    # time of infection for connection is recorded for NBDA
        }
      }
    }
    inf.status <- cbind(inf.status, c(tmp.inf.status[1:N]))    # bind infection statuses from this time step to the infection status matrix
    dat <- rbind(dat, c(tmp.dat[1:4]))    # bind SIR status to SIR dataframe
    time = time + 1    # update time step
  }
  
## model code for running non-social transmission model
} else if(model.type == 0) {
  ## run model until no more infectious individuals are present in population (0mod5==0, 5mod5==0, 2mod5!=3)
  while(sum(inf.status[,(time+1)] %% 5) > 0)
  {
    tmp.inf.status <- inf.status[,(time+1)]    # make a temporary copy of the state of the population at previous time step for updating new infection statuses
    tmp.dat <- dat[(time+1),]; tmp.dat[,1] <- time     # make a temporary copy of a dataframe row for SIR model updating
    for(i in 1:N)    # for all individuals in population
    {
      ## if individual is infected and has reached the end of the infectious period...
      if(!is.na(infection.times[i]) && (infection.times[i] + Inf.period) == (time+1))
      {
        tmp.inf.status[i] <- 5    # individual becomes recovered
        tmp.dat[,3] <- tmp.dat[,3] - 1; tmp.dat[,4] <- tmp.dat[,4] + 1    # numbers of I & R updated
      ## if individual is not infected and is randomly chosen to become infected (non-socially)...
      } else if(is.na(infection.times[i]) && rbinom(1, 1, p.inf) == 1)
      {
        tmp.inf.status[i] <- 2    # individual becomes infected
        tmp.dat[,2] <- tmp.dat[,2] - 1; tmp.dat[,3] <- tmp.dat[,3] + 1    # numbers of S & I updated
        infection.times[i] <- time + 1    # time of infection is recorded for NBDA
      }
    }
    inf.status <- cbind(inf.status, c(tmp.inf.status[1:N]))    # bind infection statuses from this time step to the infection status matrix
    dat <- rbind(dat, c(tmp.dat[1:4]))    # bind SIR status to SIR dataframe
    time = time + 1    # update time step
  }
}


#-------------------#
#   Graphics Code   #
#-------------------#

## re-arrange infection statuses so that the vector is read properly when displaying network graphics
inf.status <- rbind(inf.status[1,1:ncol(inf.status)],inf.status[10:19,1:ncol(inf.status)],inf.status[2,1:ncol(inf.status)],inf.status[20:29,1:ncol(inf.status)],inf.status[3,1:ncol(inf.status)],inf.status[30:38,1:ncol(inf.status)],inf.status[4:9,1:ncol(inf.status)])
inf.status <- replace(inf.status, inf.status == 5, 4)

## code for importing an edgelist of ties in the HEB 1333 class and converting this to a statnet network object
el <- read.csv("undirHEBV.csv",header=FALSE)
el[,1] <- as.character(el[,1])
el[,2] <- as.character(el[,2])
n <- network(el,matrix.type="edgelist",directed=FALSE)

## remember coordinates for the network graphic, otherwise this is recalculated every iteration and makes the graph impossible to follow
gcoord <- gplot.layout.kamadakawai(n, layout.par=NULL)    

## set graphing window size, with plotting corrections for Windows OS [Adapted from code by N Hardt]
if(.Platform$OS.type=="windows") {
	windows(width=12, height=6)
} else {
	quartz("Simulated Transmission of HEBV through HEB 1333", width=12, height=6)
}

for (i in 1:ncol(inf.status))    # for all time steps
{
  par(mfrow=c(1,2)); par(new=T)    # set up two side-by-side graphing boxes
  
  ## plot the network graphic on the left-hand side, with colors of nodes corresponding to infection status
  gplot(n, gmode="graph", coord=gcoord, vertex.col=inf.status[,i], cex=4, edge.col=1, displaylabels=FALSE)
  
  ## plot SIR over time on the right-hand side [Adapted from Bellan et al. 2012]
  plot(dat$time[1:i], dat$inf[1:i], type="s", col="red", ylim = c(0,45), xlim = c(0,ncol(inf.status)), lwd = 1, bty="n", lty = 1, xlab = "Time Since Outbreak (hours)", ylab = "# People")
  axis(2, at = seq(0,60, by = 20))
  lines(dat$time[1:i], dat$sus[1:i], type = "s", col = "black", lwd = 1, lty = 1)
  lines(dat$time[1:i], dat$rec[1:i], type = "s", col = "blue", lwd = 1, lty = 1)
  legend("topright", c("Susceptible", "Infected", "Recovered"), col=c("black","red","blue"), lwd = 1, cex=0.8)
  
  Sys.sleep(.1)    # wait 0.1 seconds before updating graphs
}


#-------------------#
#    Export Code    #
#-------------------#

## set up infection times into format that can be read by NBDA and export with user defined title
infection.times <- t(infection.times)
fp <- UserFilenameInput()
write.table(infection.times, fp, sep=",", col.names=FALSE, row.names=FALSE)