#' LGC: Need to add missing value functionality !! Needs a code check!
#' @keywords internal
stable_pcp <- function(D, lambda, mu, LOD = -Inf, non_negative = T, max_iter = 10000, verbose = F) {

  # 1. Initialize variables:
  # 1a. Matrix dimensions:
  n <- nrow(D)
  p <- ncol(D)

  # 1b. Process LOD input (convert LOD vector to matrix if needed, so it multiplies correctly):
  if (any(class(LOD) == "list")) LOD <- unlist(LOD)
  if (is.vector(LOD) && length(LOD) != 1) LOD <- matrix(LOD, n, p, T)

  # 1c. Omega = mask of observed values, for handling missingness later:
  Omega <- !is.na(D)
  D[!Omega] <- 0

  # 1d. L matrix primal variables (L3 carries non-negativity constraint in primal):
  L1 <- matrix(0, n, p)
  L2 <- matrix(0, n, p)
  L3 <- matrix(0, n, p)

  # 1e. S matrix primal variables:
  S1 <- matrix(0, n, p)
  S2 <- matrix(0, n, p)

  # 1f. Dual variables (noise matrix; Z2 carries non-negativity constraint in dual):
  Z1 <- matrix(0, n, p)
  Z2 <- matrix(0, n, p)
  Z3 <- matrix(0, n, p)

  # 1g. Hard-coded optimization parameters:
  eps_abs <- 1e-04 # Epsilon absolute, for stopping criteria
  eps_rel <- 1e-03 # Epsilon relative, for stopping criteria
  rho <- 0.1 # Augmented Lagrangian coefficient (rate)

  # 1h. Flag identifying convergence; vector of objective values for user
  converged <- F
  objective <- vector("numeric", max_iter)

  # 2. ADMM-splitting iterations:
  for (i in 1:max_iter) {

    # 2a. Update 1st primal variables (L1, S1):
    nuc <- prox_nuclear((L2 + L3 - (Z1 + Z2)/rho)/2, 1/2/rho)
    L1 <- nuc[[1]]
    L1_nuclear_norm <- nuc[[2]]
    S1 <- prox_l1(S2 - Z3/rho, lambda/rho)

    # 2b. Update 2nd primal variables (L2, L3, S2):
    L2_old <- L2
    S2_old <- S2
    L3_old <- L3

    L2_opt1 <- (mu*rho*D + (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1)/(2*mu*rho + rho^2)
    L2_opt2 <- L1 + Z1/rho
    L2_opt3 <- (mu*rho*LOD + (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1)/(2*mu*rho + rho^2)
    L2_opt4 <- ((mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1)/(2*mu*rho + rho^2)

    L2_new <- L2_opt1*(D >= 0) +
      L2_opt2*((D < 0) & ((L2 + S2) >= 0) & ((L2 + S2) <= LOD)) +
      L2_opt3*((D < 0) & ((L2 + S2) > LOD)) +
      L2_opt4*((D < 0) & ((L2 + S2) < 0))

    S2_opt1 <- (mu*rho*D + (mu + rho)*Z3 - mu*Z1 + (mu + rho)*rho*S1 - mu*rho*L1)/(2*mu*rho + rho^2)
    S2_opt2 <- S1 + Z3/rho
    S2_opt3 <- (mu*rho*LOD + (mu + rho)*Z3 - mu*Z1 + (mu + rho)*rho*S1 - mu*rho*L1)/(2*mu*rho + rho^2)
    S2_opt4 <- ((mu + rho)*Z3 - mu*Z1 + (mu + rho)*rho*S1 - mu*rho*L1)/(2*mu*rho + rho^2)

    S2 <- S2_opt1*(D >= 0) +
      S2_opt2*((D < 0) & ((L2 + S2) >= 0) & ((L2 + S2) <= LOD)) +
      S2_opt3*((D < 0) & ((L2 + S2) > LOD)) +
      S2_opt4*((D < 0) & ((L2 + S2) < 0))

    L2 <- L2_new
    if (non_negative) L3 <- pmax(L1 + Z2/rho, 0, na.rm = T)

    # 2c. Update dual variables (Z1, Z2, Z3):
    Z1 <- Z1 + rho*(L1 - L2)
    if (non_negative) Z2 <- Z2 + rho*(L1 - L3)
    Z3 <- Z3 + rho*(S1 - S2)

    # 2d. Calculate primal & dual residuals:
    res_primal <- sqrt(norm(L1 - L2, "F")^2 + non_negative*norm(L1 - L3, "F")^2 + norm(S1 - S2, "F")^2)
    res_dual <- rho*sqrt(norm(L2 + L3 - L2_old - L3_old, "F")^2 + norm(S2 - S2_old, "F")^2)

    # 2e. Update rho:
    if (res_primal > 10*res_dual) {
      rho <- 2*rho
    } else if (res_dual > 10*res_primal) {
      rho <- rho/2
    }

    # 2f. Calculate objective:
    objective[i] <- L1_nuclear_norm + lambda*sum(abs(S1)) + mu*loss_lod(D, L2 + S2, LOD)
    if (verbose) cat(paste0("Iteration: ", i, " Obj: ", round(objective[i], 5)))

    # 2g. Update & check stopping criteria:
    thresh_primal <- eps_abs*sqrt((4 + non_negative)*n*p) + eps_rel*max(sqrt((1 + non_negative)*norm(L1, "F")^2 + norm(S1, "F")^2), sqrt(norm(L2, "F")^2 + norm(L3, "F")^2 + norm(S2, "F")^2))
    thresh_dual <- eps_abs*sqrt((2 + non_negative)*n*p) + eps_rel*sqrt(norm(Z1 + Z2, "F")^2 + norm(Z3, "F")^2)
    if (res_primal < thresh_primal && res_dual < thresh_dual) {
      converged <- T
      if (verbose) cat(paste0("Converged in ", i, " iterations."))
      break
    }

  }

  # 3. Wrap up & return:
  if (!converged && verbose) warning(paste0("Maximum iterations reached (", max_iter, "). PCP did not converge."))

  L_final <- (L1 + L2 + L3)/(2 + non_negative)
  S_final <- (S1 + S2)/2
  if (non_negative) L_final[L_final < 0] <- 0

  list(L = L_final, S = S_final, final_iter = i, objective = objective, converged = converged)

}
