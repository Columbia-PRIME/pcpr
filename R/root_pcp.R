#' Needs a code check!
#' @keywords internal
root_pcp <- function(D, lambda, mu, LOD = -Inf, non_negative = T, max_iter = 10000, verbose = F) {

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

  # 1d. mask_above/below_LOD = entries that are observed and above/below LOD:
  mask_above_LOD <- Omega & (D >= LOD)
  mask_below_LOD <- Omega & (D < LOD)

  # 1e. Hard-coded optimization parameters:
  eps_abs <- 1e-06 # Epsilon absolute, for stopping criteria
  eps_rel <- 1e-06 # Epsilon relative, for stopping criteria
  rho <- 0.1 # Augmented Lagrangian coefficient (rate)

  # 1f. L matrix primal variables (L3 carries non-negativity constraint in primal):
  L1 <- matrix(0, n, p)
  L2 <- matrix(0, n, p)
  L3 <- matrix(0, n, p)

  # 1g. S matrix primal variables:
  S1 <- matrix(0, n, p)
  S2 <- matrix(0, n, p)

  # 1h. Noise/error matrix primal variable:
  Z  <- matrix(0, n, p)

  # 1i. Dual variables (noise matrix; Y4 carries non-negativity constraint in dual):
  Y1 <- matrix(0, n, p)
  Y2 <- matrix(0, n, p)
  Y3 <- matrix(0, n, p)
  Y4 <- matrix(0, n, p)

  # 1j. Flag identifying convergence; vector of objective values for user
  converged <- F
  objective <- vector("numeric", max_iter)

  # 2. ADMM-splitting iterations:
  for (i in 1:max_iter) {

    # 2a. Update 1st primal variables (L1, S1):
    nuc <- prox_nuclear((L2 + L3 - (Y1 + Y4)/rho)/2, 1/2/rho)
    L1 <- nuc[[1]]
    L1_nuclear_norm <- nuc[[2]]
    S1 <- prox_l1(S2 - Y2/rho, lambda/rho)

    # 2b. Update 2nd primal variables (L2, L3, S2):
    L2_old <- L2
    S2_old <- S2
    L3_old <- L3

    temp <- L2 + S2 - Y3/rho
    temp_D <- D*mask_above_LOD + temp*(mask_below_LOD & (temp >= 0) & (temp <= LOD)) + ifelse(mask_below_LOD & (temp > LOD), LOD, 0) #LOD*(mask_below_LOD & (temp > LOD))
    Z <- prox_frobenius(temp - temp_D, mu/rho) + temp_D
    if (non_negative) L3 <- pmax(L1 + Y4/rho, 0)

    L2_obs <- Omega/3*(2*L1 - S1 + Z + (2*Y1 - Y2 + Y3)/rho)
    L2_unobs <- (1 - Omega)*(L1 + Y1/rho)
    L2 <- L2_obs + L2_unobs

    S2_obs <- Omega/3*(2*S1 - L1 + Z + (2*Y2 - Y1 + Y3)/rho)
    S2_unobs <- (1 - Omega)*(S1 + Y2/rho)
    S2 <- S2_obs + S2_unobs

    # 2c. Update dual variables (Y1, Y2, Y3, Y4):
    Y1 <- Y1 + rho*(L1 - L2)
    Y2 <- Y2 + rho*(S1 - S2)
    Y3 <- Y3 + rho*Omega*(Z - (L2 + S2))
    if (non_negative) Y4 <- Y4 + rho*(L1 - L3)

    # 2d. Calculate primal & dual residuals:
    res_primal <- sqrt(norm(L1 - L2, "F")^2 + non_negative*norm(L1 - L3, "F")^2 + norm(S1 - S2, "F")^2 + norm(Omega*(Z - L2 - S2), "F")^2)
    res_dual <- rho*sqrt(norm(L2 + L3 - L2_old - L3_old, "F")^2 + norm(S2 - S2_old, "F")^2 + norm(Omega*(L2 - L2_old + S2 - S2_old), "F")^2)

    # 2e. Update rho:
    if (res_primal > 10*res_dual) {
      rho <- 2*rho
    } else if (res_dual > 10*res_primal) {
      rho <- rho/2
    }

    # 2f. Calculate objective:
    objective[i] <- L1_nuclear_norm + lambda*sum(abs(S1)) + mu*sqrt(loss_lod(D, L2 + S2, LOD))
    if (verbose) cat(paste0("Iteration: ", i, " Obj: ", round(objective[i], 5)))

    # 2g. Update & check stopping criteria:
    thresh_primal <- eps_abs*sqrt((5 + non_negative)*n*p) + eps_rel*max(sqrt((1 + non_negative)*norm(L1, "F")^2 + norm(S1, "F")^2 + norm(Z, "F")^2), sqrt(norm(L2, "F")^2 + norm(L3, "F")^2 + norm(S2, "F")^2 + norm(Omega*(L2 + S2), "F")^2))
    thresh_dual <- eps_abs*sqrt((3 + non_negative)*n*p) + eps_rel*sqrt(norm(Y1 + Y4, "F")^2 + norm(Y2, "F")^2 + norm(Y3, "F")^2)
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
