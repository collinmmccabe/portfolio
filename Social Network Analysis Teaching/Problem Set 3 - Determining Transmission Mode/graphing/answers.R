setwd("./")

epidat <- read.csv("epidemiology_data.csv")

N.susceptibles <- nrow(epidat); N.susceptibles

inf.period <- epidat$infection_end - epidat$infection_begin
avg.inf.period <- mean(inf.period); avg.inf.period

cont <- read.csv("contact_data.csv")
cont.dist <- table(factor(cont$source_id, levels = epidat$id))
Rnought <- mean(cont.dist); Rnought