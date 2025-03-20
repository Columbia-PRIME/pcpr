#' Project matrix to rank r
#'
#' @description
#' `proj_rank_r()` implements a best (i.e. closest) rank-`r` approximation of
#'  an input matrix.
#'
#' This is computed via a simple truncated singular value decomposition (SVD),
#' retaining the first `r` leading singular values/vectors of `D`.
#' This is equivalent to solving the following optimization problem:
#' \eqn{min ||X-D||_F s.t. rank(X) <= r}, where `X` is the approximated solution
#' and `D` is the input matrix.
#'
#' `proj_rank_r()` is used to iteratively model the low-rank `L` matrix in the
#' non-convex PCP function [rrmc()], providing a non-convex replacement for the
#' [prox_nuclear()] method used in the convex PCP function [root_pcp()].
#'
#' Intuitively, `proj_rank_r()` can also be thought of as providing a PCA
#' estimate of a rank-`r` matrix `L` from observed data `D`.
#'
#' @param D The input data matrix.
#' @param r The rank that `D` should be projected/truncated to.
#'
#' @returns The best rank-`r` approximation to `D` via a truncated SVD.
#'
#' @seealso [prox_nuclear()]
#' @examples
#' # Simulating a simple dataset D with the sim_data() function.
#' # The dataset will be a 10x5 matrix comprised of:
#' # 1. A rank-1 component as the ground truth L matrix; and
#' # 2. A dense Gaussian noise component corrupting L, making L full-rank
#' data <- sim_data(10, 5, 1, numeric(), 0.01)
#' # The observed matrix D is full-rank, while L is rank-1:
#' data.frame("D_rank" = matrix_rank(data$D), "L_rank" = matrix_rank(data$L))
#' before_proj_err <- norm(data$D - data$L, "F") / norm(data$L, "F")
#' # Projecting D onto the nearest rank-1 approximation, X, via proj_rank_r()
#' X <- proj_rank_r(data$D, r = 1)
#' after_proj_err <- norm(X - data$L, "F") / norm(data$L, "F")
#' proj_v_obs_err <- norm(X - data$D, "F") / norm(data$D, "F")
#' data.frame(
#'   "Observed_error" = before_proj_err,
#'   "Projected_error" = after_proj_err,
#'   "Projected_vs_observed_error" = proj_v_obs_err
#' )
#' @export
proj_rank_r <- function(D, r) {
  # Checking simple cases
  if (ncol(D) < r || nrow(D) < r) stop("r > matrix rank")
  if (ncol(D) == r || nrow(D) == r) return(D)
  # Singular value decomposition of D
  USV <- svd(D)
  # Truncate singular values yielding truncated diagonal matrix Sigma
  singular_values <- USV$d
  singular_values[(r + 1):min(dim(D))] <- 0
  Sigma_trunc <- diag(singular_values)
  # Get rank-r approximation to D using Sigma_trunc and return
  X <- USV$u %*% Sigma_trunc %*% t(USV$v)
  X
}
