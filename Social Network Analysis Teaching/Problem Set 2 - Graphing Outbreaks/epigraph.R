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

# Plot the Graph
plot(dat$tt, dat$inf, type="s", col="red", ylim = c(0,40), xlim = c(0,220), lwd = 1, bty="n", lty = 1, xlab = "Time Since Outbreak (hours)", ylab = "# People", main="Actual Disease Spread")
axis(2, at = seq(0,60, by = 20))
lines(dat$tt, dat$sus, type = "s", col = "black", lwd = 1, lty = 2)
lines(dat$tt, dat$rec, type = "s", col = "blue", lwd = 1, lty = 3)