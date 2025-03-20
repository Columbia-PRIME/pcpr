#' Proximal gradient method for the Frobenius norm
#'
#' `prox_frobenius()` implements the proximal gradient method for the Frobenius
#' norm. The proximal gradient method is used to solve non-differentiable convex
#' optimization problems. In [root_pcp()], this thresholding minimizes the
#' square root of the sum of squared error, where the error is defined as
#' `Z = D - L - S`. **This is an internal function needed by [root_pcp()]. It is
#' not expected that users should require access to this function.**
#'
#' @param Z The input error/noise matrix, `Z = D - L - S`.
#' @param c The amount by which the prox Frobenius method penalizes `Z`.
#'
#' @returns The thresholded error/noise matrix.
#'
#' @seealso [prox_l1()], [prox_nuclear()]
#' @keywords internal
prox_frobenius <- function(Z, c) {
  n <- norm(Z, "F")
  if (n <= c) {
    X <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  } else {
    X <- (1 - c / n) * Z
  }
  X
}

#' Proximal gradient method for the L1 norm
#'
#' `prox_l1()` implements the proximal gradient method for the L1 norm.
#' In [root_pcp()], this soft thresholding encourages the `S` matrix to be
#' sparse. The proximal gradient method is used to solve non-differentiable
#' convex optimization problems. **This is an internal function needed by
#' [root_pcp()]. It is not expected that users should require access to this
#' function.**
#'
#' @param S The input sparse matrix.
#' @param c The amount by which the prox L1 method penalizes `S`.
#'
#' @returns The thresholded sparse matrix.
#'
#' @seealso [hard_threshold()], [prox_frobenius()], [prox_nuclear()]
#' @keywords internal
prox_l1 <- function(S, c) {
  sign(S) * pmax(abs(S) - c, 0)
}

#' Proximal gradient method for the nuclear norm
#'
#' `prox_nuclear()` implements the proximal gradient method for the nuclear
#' norm. The nuclear norm is equivalent to the L2 norm of the singular values of
#' the matrix. This singular value thresholding encourages the low-rank `L`
#' matrix to be low-rank in [root_pcp()]. The proximal gradient method is used
#' to solve non-differentiable convex optimization problems. **This is an
#' internal function needed by [root_pcp()]. It is not expected that users
#' should require access to this function.**
#'
#' @param L The input low-rank matrix.
#' @param c The amount by which the prox nuclear method penalizes `L`.
#'
#' @returns A list containing:
#' * `X`: The thresholded low-rank matrix.
#' * `X_nuclear_norm`: The sum of the absolute values of the thresholded
#'   singular values (used in objective function).
#'
#' @seealso [proj_rank_r()], [prox_frobenius()], [prox_l1()]
#' @keywords internal
prox_nuclear <- function(L, c) {
  USV <- svd(L)
  s <- sign(USV$d) * pmax(abs(USV$d) - c, 0)
  X <- USV$u %*% diag(s) %*% t(USV$v)
  X_nuclear_norm <- sum(abs(s))
  list(X = X, X_nuclear_norm = X_nuclear_norm)
}
