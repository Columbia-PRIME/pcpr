#' Prox nuclear
#'
#' \code{prox_nuclear} implements the proximal gradient method for the nuclear norm.
#' The nuclear norm is equivalent to the L2 norm of the singular values of the matrix.
#' This singular value thresholding encourages the \code{L} matrix to be low rank.
#' The proximal gradient method is used to solve non-differentiable convex optimization problems.
#' This is used in \code{stablePCP} and \code{rootPCP}.
#'
#' @param Y The \code{L} matrix.
#' @param c The amount by which the prox nuclear method penalizes \code{Y}.
#'
#' @return \code{prox_nuclear} returns two objects.
#' \code{X} is the thresholded \code{L} matrix.
#' \code{nuclearX} is the sum of the absolute values of the thresholded singular values.
#' This goes into the objective function.
prox_nuclear <- function(Y,c) {

  USV <- svd(Y)
  U <- USV$u
  S <- USV$d
  V <- USV$v

  S_new <- sign(S) * pmax(abs(S) - c, 0)
  X <- U %*% diag(S_new) %*% t(V)
  nuclearX  <- sum(abs(S_new))

  return(list(X = X, nuclearX = nuclearX))
}

#' Prox L1
#'
#' \code{prox_l1} implements the proximal gradient method for the L1 norm.
#' This soft thresholding encourages the \code{S} matrix to be sparse.
#' The proximal gradient method is used to solve non-differentiable convex optimization problems.
#' This is used in \code{stablePCP} and \code{rootPCP}.
#'
#' @param Y The \code{S} matrix.
#' @param c The amount by which the prox L1 method penalizes \code{Y}.
#'
#' @return The thresholded \code{S} matrix.
#'
prox_l1 <- function(Y, c) {
  X <- sign(Y) * pmax(abs(Y) - c, 0)
  return(X)
}

#' Prox Frobenius
#'
#' \code{prox_fro} implements the proximal gradient method for the Frobenius norm.
#' This thresholding minimizes the square root of the sum of squared error.
#' The proximal gradient method is used to solve non-differentiable convex optimization problems.
#' This is only used in \code{rootPCP} because the squared Frobenius error in \code{stablePCP} is differentiable.
#'
#' @param X The error matrix, \code{D-L-S}
#' @param c The amount by which the prox Frobenius method penalizes \code{X}.
#'
#' @return The thresholded error matrix.
#'
prox_fro <- function(X,c) {

  n = norm(X,'F')

  if (n <= c) {Y = matrix(0, nrow = nrow(X), ncol = ncol(X))}
  else {Y = (1 - c/n) * X}

  return(Y)
}

#' Loss for stablePCP-LOD
#'
#' \code{loss_lod} includes the LOD-specific penalty terms to compute the loss on the squared error term \code{L+S-D}.
#'
#' @param X The predicted value of \code{L + S} at the current iteration.
#' @param D The original dataset.
#' @param LOD The LOD. May be a scalar, vector (\code{length(LOD) = ncol(D)}), or matrix (\code{dim(LOD) == dim(D)}).
#'
#' @return Scalar value used to calculate loss in objective function.
#'
loss_lod <- function(X, D, LOD) {
  X_lod <- (X - D)   * (D >= 0) +
    (X - LOD)  * ((D < 0) & (X > LOD)) +
    X        * ((D < 0) & (X < 0))

  sum(X_lod^2) / 2
}

#' Project rank
#' Non-convex replacement for nuclear norm
#'
#' @keywords internal
old_proj_rank_r = function(Y, r) {

  if (ncol(Y) < r) stop("r > matrix rank")
  if (ncol(Y) == r) {return(Y)}

  USV <- svd(Y)
  U <- USV$u
  S <- USV$d
  V <- USV$v

  s = S
  
  # used to read:
  #if (length(S)-1 == r) {
  #  s[r]  = 0
  #} else {
  #  s[(r+1):length(s)] = 0
  #}
  # now reads:
  s[(r+1):length(s)] = 0

  S_new  = diag(s)
  X = U %*% S_new %*% t(V)
  return(X)
}

#' @export
proj_r_partial <- function(Y, m1, m2, r, numIter = 10, X11_L = NULL, verbose = FALSE) {
  #  % solve the optimization problem 
  #%
  #%   min_{rank(X) <= r} .5 || P_Gamma^c [Y - X] ||_F^2
  #%
  #% via the iteration 
  #% 
  #%   X^+ = P_r  [  X11  Y12  ]
  #%              [  Y21  Y22  ]
  #% 
  #% 
  #% This method assumes that entries missing are the upper left block 
  #%    { 1, ..., m1 } x { 1, ..., m2 }
  #%  
  #%  

  #% parameters
  Y[is.na(Y)] <- 0

  M1 <- nrow(Y)
  M2 <- ncol(Y)

  #% initialize
  Y12 <- Y[1:m1, (m2+1):M2]
  Y22 <- Y[(m1+1):M1, (m2+1):M2]
  Y21 <- Y[(m1+1):M1, 1:m2]

  UsV <- svd(Y22)
  U <- UsV$u
  s <- UsV$d
  V <- UsV$v

  Y22_pinv <- V[, 1:r] %*% solve(diag(s)[1:r, 1:r], tol = 1e-18) %*% t(U[, 1:r])

  X11 <- Y12 %*% Y22_pinv %*% Y21

  X <- Y

  X[1:m1, 1:m2] <- X11

  X <- proj_rank_r(X, r)
  
  #X11_init <- X11

  #for (i in 1:numIter) {
  #  H <- Y
  #  H[1:m1, 1:m2] <- X11

  #  X <- proj_rank_r(H, r)

  #  X11_new <- X[1:m1, 1:m2]

  #  if (is.null(X11_L)) {
  #    step <- norm(X11_new - X11, "F")
  #    if (verbose) cat(paste("Iter:", i, "Step: ", step, "\n"))
  #  } else {
  #    step <- norm(X11_new - X11_L)
  #    if (verbose) cat(paste("Iter:", i, "Loss: ", step, "\n"))
  #  }
    
  #  X11 <- X11_new

  #}

  return(X)
}

#' HT
#' Hard-thresholding for the sparse matrix
#' @keywords internal
HT = function(Y, c) {
  X = Y
  X[abs(X) < c] = 0
  return(X)
}