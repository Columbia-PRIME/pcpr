#' Prox Frobenius
#'
#' \code{prox_frobenius} implements the proximal gradient method for the Frobenius norm.
#' This thresholding minimizes the square root of the sum of squared error.
#' The proximal gradient method is used to solve non-differentiable convex optimization problems.
#' This is only used in \code{root_pcp}, since the squared Frobenius error in \code{stable_pcp} is differentiable.
#'
#' @param Z The input error/noise matrix, \code{Z = D - L - S}.
#' @param c The amount by which the prox Frobenius method penalizes \code{Z}.
#'
#' @return The thresholded error/noise matrix.
#'
#' @seealso \code{\link{prox_l1}}, \code{\link{prox_nuclear}}
#' @examples
#' @keywords internal
prox_frobenius <- function(Z, c) {

  n <- norm(Z, "F")

  if (n <= c) {
    X <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  } else {
    X <- (1 - c/n) * Z
  }

  X

}

#' Prox L1
#'
#' \code{prox_l1} implements the proximal gradient method for the L1 norm.
#' This soft thresholding encourages the \code{S} matrix to be sparse.
#' The proximal gradient method is used to solve non-differentiable convex optimization problems.
#' This is used in \code{stable_pcp} and \code{root_pcp}.
#'
#' @param S The input sparse matrix.
#' @param c The amount by which the prox L1 method penalizes \code{S}.
#'
#' @return The thresholded sparse matrix.
#'
#' @seealso \code{\link{hard_thresholding}}, \code{\link{prox_frobenius}}, \code{\link{prox_nuclear}}
#' @examples
#' @keywords internal
prox_l1 <- function(S, c) {

  X <- sign(S) * pmax(abs(S) - c, 0)
  X

}

#' Prox Nuclear
#'
#' \code{prox_nuclear} implements the proximal gradient method for the nuclear norm.
#' The nuclear norm is equivalent to the L2 norm of the singular values of the matrix.
#' This singular value thresholding encourages the low-rank \code{L} matrix to be low rank.
#' The proximal gradient method is used to solve non-differentiable convex optimization problems.
#' This is used in \code{stable_pcp} and \code{root_pcp}.
#'
#' @param L The input low-rank matrix.
#' @param c The amount by which the prox nuclear method penalizes \code{L}.
#'
#' @return \code{prox_nuclear} returns a list containing the following two objects:
#' \describe{
#'    \item{\code{X}}{The thresholded low-rank matrix.}
#'    \item{\code{X_nuclear_norm}}{The sum of the absolute values of the thresholded singular values (used in objective function).}
#' }
#'
#' @seealso \code{\link{proj_rank_r}}, \code{\link{prox_frobenius}}, \code{\link{prox_l1}}
#' @examples
#' @keywords internal
prox_nuclear <- function(L, c) {

  USV <- svd(L)

  s <- sign(USV$d) * pmax(abs(USV$d) - c, 0)
  X <- USV$u %*% diag(s) %*% t(USV$v)
  X_nuclear_norm <- sum(abs(s))

  list(X = X, X_nuclear_norm = X_nuclear_norm)

}

#' Loss function for PCP's limit of detection (LOD) penalty.
#'
#' \code{loss_lod} includes the LOD-specific penalty terms to compute the loss on the squared error term \code{L + S - D}.
#'
#' @param D The original data matrix.
#' @param X The predicted value of \code{L + S} at the current iteration.
#' @param LOD The LOD as a matrix (\code{dim(LOD) == dim(D)}).
#'
#' @return Scalar value used to calculate loss in objective function.
#'
#' @seealso \code{\link{stable_pcp}}
#' @examples
#' @keywords internal
loss_lod <- function(D, X, LOD) {

  X_lod <- (X - D)*(D >= 0) + (X - LOD)*((D < 0) & (X > LOD)) + X*((D < 0) & (X < 0))
  sum(X_lod^2) / 2

}
