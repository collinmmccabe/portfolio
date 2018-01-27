dat$susprop <-  dat$sus / rowSums(dat[,2:4])
var2_est <- var1_est*dat$susprop