# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

############################################################
## Original PCP function
############################################################

soft_thresholding <- function(v, lambda) {
  myzero <- matrix(data = 0, ncol = ncol(v), nrow = nrow(v))
  w <- sign(v) * pmax(abs(v) - lambda, myzero)
  w
}

############################################################

soft_thresholding_diag <- function(v, lambda) {
  myzero <- vector("numeric", length = length(v))
  w <- sign(v) * pmax(abs(v) - lambda, myzero)
  w
}

############################################################

singular_value_threshold <- function(M, lambda) {

  USV <- svd(M)
  U <- USV$u
  S <- USV$d
  V <- USV$v

  N <- U %*% diag(soft_thresholding_diag(S, lambda)) %*% t(V)

  v  <- sum(soft_thresholding_diag(S, lambda))

  svt <- list(N = N, v = v)
  svt
}

############################################################

original_pcp <- function(D, lambda, mu) {

  D <- as.matrix(D)
  m <- nrow(D)
  n <- ncol(D)

  S <- matrix(0, nrow = m, ncol = n)
  L <- matrix(0, nrow = m, ncol = n)

  iter <- 0
  MAX_ITER <- 5000
  done <- FALSE

  # Convergence Thresholds
  LOSS_THRESH <- 1e-4
  loss <- vector("numeric", MAX_ITER)

  while (!done) {

    iter <- iter + 1

    svt <- singular_value_threshold((D - S), 1/mu)
    L <- svt[[1]] #svt$N
    v <- svt[[2]]

    S <- soft_thresholding((D - L), lambda/mu)

    obj <- v + lambda * sum(abs(S)) + (mu/2) * norm((D - L - S), type = "F")^2
    loss[iter] <- obj

    print(pate0(iter, " Obj: ", obj))

    if (iter >= MAX_ITER |
        (iter != 1) && (abs(loss[iter-1] - loss[iter]) < LOSS_THRESH)) {done <- TRUE}

  }
  list(L = L, S = S)
}

############################################################
## PCP-LOD, new function with penalty for values <LOD
############################################################

# Prox L1 norm function, soft thresholding
# if Y < c (threshold), push to zero
prox_l1 <- function(Y, c) {

  myzero <- matrix(data = 0, ncol = ncol(Y), nrow = nrow(Y))
  X <- sign(Y) * pmax(abs(Y) - c, myzero, na.rm = TRUE)
  X
}

############################################################

# Prox nuclear norm function, L1 norm of the singular values
# This encourages matrix to be low rank by pushing SV to zero (sparse)
prox_nuclear <- function(Y,c) {

  USV <- svd(Y)
  U <- USV$u
  S <- USV$d
  V <- USV$v

  myzero <- vector("numeric", length = length(S))
  S_new <- sign(S) * pmax(abs(S) - c, myzero, na.rm = TRUE)
  # Threshold the singular values, if SV < c, push it to zero

  X <- U %*% diag(S_new) %*% t(V)
  # this t(.) is for all scalar, vector, matrices
  # % X is the truncation of the original
  # % Multiply the thresholded SVD components back together

  nuclearX  <- sum(abs(S_new))
  # This is the L1 norm of the truncated singular values
  # Goes into the loss function

  list(X = X, nuclearX = nuclearX)
}

############################################################

# is same function for convergence criteria
# is the difference among matrices > noise threshold?
## if TRUE, keep iterating, if FALSE, end

# Compares L1, L2, L3 OR S1, S2
# Need ... for function to handle different number of inputs
# length(varargin) gives the number of function input arguments given in the call
# for L1, L2, L3, THREE comparisons, L1/L2, L1/L3, and L2/L3
# for S1, S2, one comparison
is_same <- function(SAME_THRESH, ...) {
  flag <- TRUE
  varargin <- list(...)
  if (length(varargin) == 2) {
    if (max(abs(varargin[[1]] - varargin[[2]])) > SAME_THRESH) {
      flag <- FALSE
    }
  }
  else if (length(varargin) == 3) {
    if ((max(abs(varargin[[1]] - varargin[[2]])) > SAME_THRESH) |
        (max(abs(varargin[[1]] - varargin[[3]])) > SAME_THRESH) |
        (max(abs(varargin[[2]] - varargin[[3]])) > SAME_THRESH)) {
      flag <- FALSE
    }
  }
  flag
}

############################################################

# loss_lod function (only used in the loss function)
loss_lod <- function(X, D, LOD) {
  # % D is the original data
  # % X is the new thing (L + S)
  # # LOD is the LOD
  X_lod <- (X - D)   * (D >= 0) +
    (X - LOD)  * ((D < 0) & (X > LOD)) +
    X        * ((D < 0) & (X < 0)) #+
  #(X - D)    * (D < LOD & (X > 0 && X <= LOD )) # or should it be zero

  l <- sum(X_lod^2) / 2
  # % L2 norm

  # % Any D_ij < 0 AND X_ij < LOD AND > 0 are treated as equal

  # % Minimize discrepancy for valid data
  # % Want to shrink negative things
  l
}

############################################################

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

  m <- nrow(D)
  n <- ncol(D)
  rho <- 1 # Augmented Lagrangian coefficient (rate)

  L1 <- matrix(0, m, n)
  L2 <- matrix(0, m, n)
  L3 <- matrix(0, m, n)

  S1 <- matrix(0, m, n)
  S2 <- matrix(0, m, n)

  Z1 <- matrix(0, m, n)
  Z2 <- matrix(0, m, n)
  Z3 <- matrix(0, m, n)

  # Max iteration
  MAX_ITER <- 1000

  # Convergence Thresholds
  LOSS_THRESH <- 1e-5
  SAME_THRESH <- 1e-4

  # if (is.vector(LOD)) {
  #   tf = ifelse(D < 0, TRUE, FALSE)
  #   LOD = t(t(tf) * LOD)
  # }

  if (is.vector(LOD)) {
    empty = matrix(1, nrow = nrow(D), ncol = ncol(D))
    LOD = t(t(empty) * LOD)
  } # This converts a vector LOD to a matrix, so that it multiplies correctly

  loss <- vector("numeric", MAX_ITER)

  for (i in 1:MAX_ITER) {

    nuc <- prox_nuclear( ((L2 + L3 - (Z1 + Z2)/rho)/2), 1/2/rho)
    L1 <- nuc[[1]]
    nuclearL1 <- nuc[[2]] #nuclearX
    # % L, Z, S all start at zero, and change each iteration
    # % Prox_nuc is singular value thresholding
    # % L is low rank matrix

    S1 <- prox_l1((S2 - Z3/rho), lambda/rho)
    # % S is sparse matrix

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

    print(paste0(i, " Obj: ", loss[i]))

    if ((i != 1) &&
        (abs(loss[i-1] - loss[i]) < LOSS_THRESH) &&
        is_same(SAME_THRESH, L1, L2, L3) &&
        is_same(SAME_THRESH, S1, S2)) {
      break} # % Convergence criteria!
  }

  L <- L3 # (L1 + L2 + L3) / 3
  S <- S1 #(S1 + S2) / 2
  list(L = L, S = S)
}

############################################################
## PCP-LOD version where sparse matrix is non-negative, too
############################################################

pcp_lod_nnS <- function(D, lambda, mu, LOD) {

  m <- nrow(D)
  n <- ncol(D)
  rho <- 1 # Augmented Lagrangian coefficient (rate)

  L1 <- matrix(0, m, n)
  L2 <- matrix(0, m, n)
  L3 <- matrix(0, m, n)

  S1 <- matrix(0, m, n)
  S2 <- matrix(0, m, n)
  S3 <- matrix(0, m, n) # ADDED

  Z1 <- matrix(0, m, n)
  Z2 <- matrix(0, m, n)
  Z3 <- matrix(0, m, n)
  Z4 <- matrix(0, m, n) # ADDED

  # Max iteration
  MAX_ITER <- 5000

  # Convergence Thresholds
  LOSS_THRESH <- 1e-5
  SAME_THRESH <- 1e-4

  if (is.vector(LOD)) {
    empty = matrix(1, nrow = nrow(D), ncol = ncol(D))
    LOD = t(t(empty) * LOD)
  } # This converts a vector LOD to a matrix, so that it multiplies correctly

  loss <- vector("numeric", MAX_ITER)

  for (i in 1:MAX_ITER) {

    nuc <- prox_nuclear(((L2 + L3 - (Z1 + Z2)/rho)/2), 1/2/rho)
    L1 <- nuc[[1]]
    nuclearL1 <- nuc[[2]] #nuclearX

    S1 <- prox_l1(((S2 + S3 - (Z3 + Z4)/rho)/2), lambda/rho) # ADDED
    #prox_l1(S2 - Z3/rho, lambda/rho)

    L2_opt1 <- (mu*rho*D     + (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)
    L2_opt2 <- L1 + Z1/rho
    L2_opt3 <- ((mu*rho*LOD + (((mu + rho)*Z1) - (mu*Z3) + ((mu + rho)*rho*L1) - (mu*rho*S1)))) / ((2*mu*rho) + (rho^2))
    L2_opt4 <- (               (mu + rho)*Z1 - mu*Z3 + (mu + rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2)

    L2 <- (L2_opt1 * (D >= 0)) +
      (L2_opt2 * ((D < 0) & (((L2 + S2) >= 0) & ((L2 + S2) <= LOD)))) +
      (L2_opt3 * ((D < 0) & (((L2 + S2) > LOD)))) +
      (L2_opt4 * ((D < 0) & (((L2 + S2) < 0))))

    S2_opt1 <- (mu*rho*D     + (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)
    S2_opt2 <- S1 + (Z3/rho)
    S2_opt3 <- (((mu*rho*LOD) + (((mu + rho)*Z3) - (mu*Z1) + ((mu + rho)*rho*S1) - (mu*rho*L1)))) / ((2*mu*rho) + (rho^2))
    S2_opt4 <- (               (mu + rho)*Z3 - (mu*Z1) + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2)

    S2 <- (S2_opt1 * (D >= 0)) +
      (S2_opt2 * (((D < 0) & ((L2 + S2) >= 0) & ((L2 + S2) <= LOD)))) +
      (S2_opt3 * (((D < 0) & ((L2 + S2) > LOD)))) +
      (S2_opt4 * (((D < 0) & ((L2 + S2) < 0))))

    L3 <- pmax(L1 + Z2/rho, 0, na.rm = TRUE)
    # % Non-Negativity constraint!

    ## ADDED
    S3 <- pmax(S1 + Z4/rho, 0, na.rm = TRUE)

    Z1 <- Z1 + rho*(L1 - L2)
    Z2 <- Z2 + rho*(L1 - L3)
    Z3 <- Z3 + rho*(S1 - S2)
    Z4 <- Z4 + rho*(S1 - S3) # ADDED
    # % Z accumulate differnces between L and L and between S and S

    loss[i] <- nuclearL1 +
      (lambda*sum(abs(S1))) +
      (mu*loss_lod((L2 + S2), D, LOD)) +
      sum(Z1*(L1 - L2)) +
      sum(Z2*(L1 - L3)) +
      sum(Z3*(S1 - S2)) +
      sum(Z4*(S1 - S3)) + # ADDED
      (rho/2 * (sum((L1-L2)^2) + sum((L1 - L3)^2) + sum((S1 - S2)^2)) + sum((S1 - S3)^2)) # ADDED
    # % The code block above takes LOD into account.

    print(paste0(i, " Obj: ", loss[i]))

    if ((i != 1) &&
        (abs(loss[i-1] - loss[i]) < LOSS_THRESH) &&
        is_same(SAME_THRESH, L1, L2, L3) &&
        is_same(SAME_THRESH, S1, S2, S3)) { # ADDED
      break} # % Convergence criteria!
  }

  L <- L3 #(L1 + L2 + L3) / 3
  S <- S3 #(S1 + S2 + S3) / 3
  list(L = L, S = S)
}

############################################################
############################################################





