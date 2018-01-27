cont.dist <- table(factor(cont$source_id, levels = epidat$id))
var1_est <- mean(cont.dist)