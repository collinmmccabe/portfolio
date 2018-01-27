setwd("../")


#---
# S-I-R Calculations
epidat <- read.csv("epidemiology_data.csv")
event.times <- c(epidat$infection_begin, epidat$infection_end)
event.times <- unique(event.times)
event.times <- event.times[order(event.times)]
event.times <- c(-12, event.times)

n0_susceptibles <- nrow(epidat)
dat <- data.frame(tt = event.times, sus = NA, inf = NA, rec = NA)
dat[1,2:4] <- c(n0_susceptibles, 0 , 0)

for(ii in 2:nrow(dat))
  {
    tt <- dat$tt[ii]
    ## Figure out how many new infections were during that event
    new.inf <- sum(epidat$infection_begin == tt, na.rm = T)
    ## Figure out how many new recoveries were during that event    
    new.rec <- sum(epidat$infection_end == tt, na.rm = T)
    ## Calculate new state variables
    dat$sus[ii] <- dat$sus[ii-1]-new.inf
    dat$inf[ii] <- dat$inf[ii-1]+new.inf-new.rec
    dat$rec[ii] <- dat$rec[ii-1]+new.rec
  }

dat <- rbind(dat, dat[nrow(dat),]) 
dat$tt[nrow(dat)] <- max(event.times)+12

# Plot the S-I-R Graph
plot(dat$tt, dat$inf, type="s", col="red", ylim = c(0,40), xlim = c(0,220), lwd = 1, bty="n", lty = 1, xlab = "Time Since Outbreak (hours)", ylab = "# People", main="Actual Disease Spread")
axis(2, at = seq(0,60, by = 20))
lines(dat$tt, dat$sus, type = "s", col = "black", lwd = 1, lty = 2)
lines(dat$tt, dat$rec, type = "s", col = "blue", lwd = 1, lty = 3)  


# ---
# R0 estimation
cont <- read.csv("contact_data.csv")
cont.dist <- table(factor(cont$source_id, levels = epidat$id))
var1_est <- mean(cont.dist); var1_est

# RE estimation
dat$susprop <-  dat$sus / rowSums(dat[,2:4])
var2_est <- var1_est*dat$susprop; var2_est


# Plot comparison of R values
plot(dat$tt, var2_est, type = "s", col = "purple", lty = 1, lwd = 1, ylab = expression(R[E]), ylim = c(0,4), xlim=c(0,220), xlab = "Time Since Outbreak (hours)", bty="n", main="Disease Spread Variable Estimates")
abline(h = var1_est, lty = 2)


# ---
# Vaccinations Calculation
p <- 1-(1/var1_est)
n_vacc <- p*4500000; n_vacc