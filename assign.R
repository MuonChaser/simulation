
assign_strata = function(data) {
A_j <- sample(c(rep(0, J/2), rep(1, J/2)))
data$A <- A_j[data$id]

z_prop_high <- 0.7
z_prop_low <- 0.3

data$Z <- rep(0, N)
for (j in 1:J) {
  idx_j <- which(data$id == j)
  n_j <- length(idx_j)
  B_j <- data$B[idx_j]
  if (A_j[j] == 1) {
    for (k in 1:K) {
      idx_jk <- idx_j[which(B_j == k)]
      n_jk <- length(idx_jk)
      n_treat <- round(z_prop_high * n_jk)
      treat_indices <- sample(idx_jk, n_treat)
      data$Z[treat_indices] <- 1
    }
  } else {
    for (k in 1:K) {
      idx_jk <- idx_j[which(B_j == k)]
      n_jk <- length(idx_jk) 
      n_treat <- round(z_prop_low * n_jk)
      treat_indices <- sample(idx_jk, n_treat)
      data$Z[treat_indices] <- 1
    }
  }
}


data$D <- ifelse(data$Z == 1 & data$A == 1, data$d_z1_a1,
                 ifelse(data$Z == 0 & data$A == 1, data$d_z0_a1,
                        ifelse(data$Z == 1 & data$A == 0, data$d_z1_a0,
                               data$d_z0_a0)))
data$Y <- ifelse(data$Z == 1 & data$A == 1, data$y_z1_a1,
                 ifelse(data$Z == 0 & data$A == 1, data$y_z0_a1,
                        ifelse(data$Z == 1 & data$A == 0, data$y_z1_a0,
                               data$y_z0_a0)))
return(data)
}


assign_no_strata = function(data) {
A_j <- sample(c(rep(0, J/2), rep(1, J/2)))
data$A <- A_j[data$id]

z_prop_high <- 0.7
z_prop_low <- 0.3

data$Z <- rep(0, N)
for (j in 1:J) {
  idx_j <- which(data$id == j)
  n_j <- length(idx_j)
  prop_treat <- ifelse(A_j[j] == 1, z_prop_high, z_prop_low)
  n_treat <- round(prop_treat * n_j)
  treat_indices <- sample(idx_j, n_treat)
  data$Z[treat_indices] <- 1
}


data$D <- ifelse(data$Z == 1 & data$A == 1, data$d_z1_a1,
                 ifelse(data$Z == 0 & data$A == 1, data$d_z0_a1,
                        ifelse(data$Z == 1 & data$A == 0, data$d_z1_a0,
                               data$d_z0_a0)))
data$Y <- ifelse(data$Z == 1 & data$A == 1, data$y_z1_a1,
                 ifelse(data$Z == 0 & data$A == 1, data$y_z0_a1,
                        ifelse(data$Z == 1 & data$A == 0, data$y_z1_a0,
                               data$y_z0_a0)))
return(data)
}