## verify_2sls_ratio.R
## Numerically verify Theorem 3:
##   CADE_hat(a) = DEY_hat(a) / DED_hat(a)  ==  beta^Y_2SLS_{1a}  from CADEreg_new
##
## The direct ratio estimator follows the formulas in the paper:
##   DED_jk(a)  = mean_D(Z=1, B=k, cluster=j) - mean_D(Z=0, B=k, cluster=j)
##   DEY_jk(a)  = mean_Y(Z=1, B=k, cluster=j) - mean_Y(Z=0, B=k, cluster=j)
##   DED_j(a)   = (1/n_j) * sum_k [ n_jk * DED_jk(a) ]
##   DEY_j(a)   = (1/n_j) * sum_k [ n_jk * DEY_jk(a) ]
##   DED(a)     = (1/J_a) * sum_{j: A_j=a} DED_j(a)
##   DEY(a)     = (1/J_a) * sum_{j: A_j=a} DEY_j(a)
##   CADE(a)    = DEY(a) / DED(a)

source("../CADEreg_new.R")
source("../assign.R")
source("../gen_data.R")

set.seed(42)

## --- 1. Generate population data ---
pop  <- gen_data(J = 10, K = 4, n_jk = 20)
data <- pop$data

## --- 2. Assign treatment for one replication ---
obs <- assign_strata(data)

## --- 3. Direct ratio estimator ---
ratio_CADE <- function(obs, a_val) {
  J_a <- sum(tapply(obs$A, obs$id, unique) == a_val)
  clusters <- levels(obs$id)

  DED_j <- numeric(length(clusters))
  DEY_j <- numeric(length(clusters))

  for (j in seq_along(clusters)) {
    idx_j <- obs$id == clusters[j]
    if (unique(obs$A[idx_j]) != a_val) next

    n_j <- sum(idx_j)
    strata <- sort(unique(obs$B[idx_j]))

    DED_jk <- numeric(length(strata))
    DEY_jk <- numeric(length(strata))
    n_jk_vec <- numeric(length(strata))

    for (s in seq_along(strata)) {
      k <- strata[s]
      idx_jk <- idx_j & obs$B == k
      n_jk_vec[s] <- sum(idx_jk)

      D1 <- mean(obs$D[idx_jk & obs$Z == 1])
      D0 <- mean(obs$D[idx_jk & obs$Z == 0])
      Y1 <- mean(obs$Y[idx_jk & obs$Z == 1])
      Y0 <- mean(obs$Y[idx_jk & obs$Z == 0])

      DED_jk[s] <- D1 - D0
      DEY_jk[s] <- Y1 - Y0
    }
    DED_j[j] <- sum(n_jk_vec * DED_jk) / n_j
    DEY_j[j] <- sum(n_jk_vec * DEY_jk) / n_j
  }

  # average over clusters with A_j == a_val
  A_j <- tapply(obs$A, obs$id, unique)
  mask <- A_j == a_val
  DED_hat <- mean(DED_j[mask])
  DEY_hat <- mean(DEY_j[mask])

  list(CADE = DEY_hat / DED_hat, DED = DED_hat, DEY = DEY_hat)
}

ratio1 <- ratio_CADE(obs, a_val = 1)
ratio0 <- ratio_CADE(obs, a_val = 0)

cat("=== Direct ratio estimator ===\n")
cat(sprintf("CADE(1): %.8f  [DEY=%.6f, DED=%.6f]\n", ratio1$CADE, ratio1$DEY, ratio1$DED))
cat(sprintf("CADE(0): %.8f  [DEY=%.6f, DED=%.6f]\n", ratio0$CADE, ratio0$DEY, ratio0$DED))

## --- 4. 2SLS estimator via CADEreg_new ---
reg_result <- CADEreg_new(obs)

cat("\n=== 2SLS estimator (CADEreg_new) ===\n")
cat(sprintf("CADE(1): %.8f\n", reg_result$CADE1))
cat(sprintf("CADE(0): %.8f\n", reg_result$CADE0))

## --- 5. Comparison ---
cat("\n=== Differences (2SLS - ratio) ===\n")
cat(sprintf("|CADE(1) diff|: %.2e\n", abs(reg_result$CADE1 - ratio1$CADE)))
cat(sprintf("|CADE(0) diff|: %.2e\n", abs(reg_result$CADE0 - ratio0$CADE)))

tol <- 1e-8
if (abs(reg_result$CADE1 - ratio1$CADE) < tol &&
    abs(reg_result$CADE0 - ratio0$CADE) < tol) {
  cat("\nVERIFICATION PASSED: 2SLS == ratio estimator (within numerical tolerance).\n")
} else {
  cat("\nVERIFICATION FAILED: differences exceed tolerance.\n")
}

## --- 6. Repeat across multiple random seeds for robustness ---
cat("\n=== Robustness check (100 replications) ===\n")
max_diff1 <- 0; max_diff0 <- 0
for (seed in 1:100) {
  set.seed(seed)
  obs_s <- assign_strata(data)
  r1 <- ratio_CADE(obs_s, 1); r0 <- ratio_CADE(obs_s, 0)
  rr <- CADEreg_new(obs_s)
  max_diff1 <- max(max_diff1, abs(rr$CADE1 - r1$CADE))
  max_diff0 <- max(max_diff0, abs(rr$CADE0 - r0$CADE))
}
cat(sprintf("Max |CADE(1) diff| over 100 seeds: %.2e\n", max_diff1))
cat(sprintf("Max |CADE(0) diff| over 100 seeds: %.2e\n", max_diff0))
if (max_diff1 < tol && max_diff0 < tol) {
  cat("All replications: VERIFICATION PASSED.\n")
} else {
  cat("Some replications: VERIFICATION FAILED.\n")
}
