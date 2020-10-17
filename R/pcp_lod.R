#' PCP-LOD function
#'
#' @export

# % If the LOD threshold LOD = 0, solve the following ADMM splitting problem:
#   % min_{L1,L2,L3,S1,S2}
# %      ||L1||_* + lambda * ||S1||_1 + mu/2 * ||L2+S2-D||_F^2 + I_{L3>=0}
# % s.t. L1 = L2
# %      L1 = L3
# %      S1 = S2.
# %
# % If LOD is not 0, replace ||L2+S2-D||_F^2 with LOD penalty.
# %
# % Below-LOD data input in D should be denoted as negative values, e.g. -1.

pcp_lod <- function(D, lambda, mu, LOD) {

  n <- nrow(D)
  p <- ncol(D)
  rho <- 0.1 # Augmented Lagrangian coefficient (rate)

  L1 <- matrix(0, n, p)
  L2 <- matrix(0, n, p)
  L3 <- matrix(0, n, p)

  S1 <- matrix(0, n, p)
  S2 <- matrix(0, n, p)

  Z1 <- matrix(0, n, p)
  Z2 <- matrix(0, n, p)
  Z3 <- matrix(0, n, p)

  # Max iteration
  MAX_ITER <- 10000

  #%%%%% BEGIN NEW %%%%%
  EPS_ABS = 1e-4
  EPS_REL = 1e-3

  flag_converge = 0
  #%%%%% END NEW %%%%%

  # Convergence Thresholds
  LOSS_THRESH <- 1e-5
  SAME_THRESH <- 1e-4

  if (is.vector(LOD)) {
    empty = matrix(1, nrow = nrow(D), ncol = ncol(D))
    LOD = t(t(empty) * LOD)
  } # This converts a vector LOD to a matrix, so that it multiplies correctly

  loss <- vector("numeric", MAX_ITER)

  #% ADMM-splitting iterations
  for (i in 1:MAX_ITER) {

    #% Update 1st primal variable (L1,S1)
    nuc <- prox_nuclear( ((L2 + L3 - (Z1 + Z2)/rho)/2), 1/2/rho)
    L1 <- nuc[[1]]
    nuclearL1 <- nuc[[2]] #nuclearX
    # % L, Z, S all start at zero, and change each iteration
    # % Prox_nuc is singular value thresholding

    S1 <- prox_l1((S2 - Z3/rho), lambda/rho)
    # % S is sparse matrix

    #% Update 2nd primal variable (L2,L3,S2)
    #%%%%% BEGIN NEW %%%%%
    L2_old = L2
    L3_old = L3
    S2_old = S2
    #%%%%% END NEW %%%%%

    L2_opt1 <- (mu*rho*D     + (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)
    L2_opt2 <- L1 + Z1/rho
    L2_opt3 <- ((mu*rho*LOD + (((mu + rho)*Z1) - (mu*Z3) + ((mu + rho)*rho*L1) - (mu*rho*S1)))) / ((2*mu*rho) + (rho^2))
    L2_opt4 <- (               (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)

    L2_new <- (L2_opt1 * (D >= 0)) +
      (L2_opt2 * (((D < 0) & ((L2 + S2) >= 0) & ((L2 + S2) <= LOD)))) +
      (L2_opt3 * (((D < 0) & ((L2 + S2) > LOD)))) +
      (L2_opt4 * (((D < 0) & ((L2 + S2) < 0))))

    S2_opt1 <- (mu*rho*D     + (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
    S2_opt2 <- S1 + (Z3/rho)
    S2_opt3 <- (((mu*rho*LOD) + (((mu + rho)*Z3) - (mu*Z1) + ((mu + rho)*rho*S1) - (mu*rho*L1)))) / ((2*mu*rho) + (rho^2))
    S2_opt4 <- (               (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)

    S2 <- (S2_opt1 * (D >= 0)) +
      (S2_opt2 * ((D < 0) & (((L2 + S2) >= 0) & ((L2 + S2) <= LOD)))) +
      (S2_opt3 * ((D < 0) & (((L2 + S2) > LOD)))) +
      (S2_opt4 * ((D < 0) & (((L2 + S2) < 0))))

    L2 <- L2_new
    L3 <- pmax(L1 + Z2/rho, 0, na.rm = TRUE)
    # % Non-Negativity constraint!

    #%%%%% BEGIN NEW %%%%%
    #% Calculate primal & dual residuals; Update rho
    res_primal = sqrt( norm(L1-L2, "F")^2 + norm(L1-L3, "F")^2 + norm(S1-S2, "F") ^2)
    res_dual = rho * sqrt( norm(L2+L3-L2_old-L3_old,'F')^2 + norm(S2-S2_old,'F')^2 )

    if (res_primal > 10 * res_dual) {
      rho = rho * 2
    } else if (res_dual > 10 * res_primal) {
        rho = rho / 2}
    #%%%%% END NEW %%%%%

    #% Update dual variable (Z1,Z2,Z3)
    Z1 <- Z1 + rho*(L1 - L2)
    Z2 <- Z2 + rho*(L1 - L3)
    Z3 <- Z3 + rho*(S1 - S2)
    # % Z accumulate differnces between L and L and between S and S

    loss[i] <- nuclearL1 +
      (lambda*sum(abs(S1))) +
      (mu*loss_lod((L2 + S2), D, LOD)) +
      sum(Z1*(L1 - L2)) +
      sum(Z2*(L1 - L3)) +
      sum(Z3*(S1 - S2)) +
      (rho/2 * (sum((L1-L2)^2) + sum((L1 - L3)^2) + sum((S1 - S2)^2)))
    # % The code block above takes LOD into account.

    #%%%%% BEGIN NEW %%%%%
    #% Check stopping criteria
    thresh_primal = EPS_ABS * sqrt(3*n*p) + EPS_REL * pmax(sqrt( norm(L1,'F')^2 * 2 + norm(S1,'F')^2 ),
                                                           sqrt( norm(L2,'F')^2 + norm(L3,'F')^2 + norm(S2,'F')^2 ))
    thresh_dual = EPS_ABS * sqrt(2*n*p) + EPS_REL * sqrt( norm(Z1+Z2,'F')^2 + norm(Z3,'F')^2 )

    if (res_primal < thresh_primal && res_dual < thresh_dual) {
      flag_converge = 1;
      print(paste0('Converged in ', i,' iterations.'))
      break}
    #%%%%% END NEW %%%%%
  }

  if( i == MAX_ITER) warning('Maximum iterations reached. PCP did not converge.')
#  print(paste0("Iteration: ", i, " Obj: ", round(loss[i], 5)))

  L <- L3 # (L1 + L2 + L3) / 3
  S <- S1 #(S1 + S2) / 2
  list(L = L, S = S)
}

############################################################
## PCP with nothing <LOD
############################################################

pcp_original <- function(D, lambda, mu) {
  pcp_lod(D, lambda, mu, LOD = 0)
  }




