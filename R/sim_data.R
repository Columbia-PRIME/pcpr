#' Simulate simple mixtures data
#'
#' @description
#' `sim_data()` generates a simulated dataset `D = L + S + Z` for
#' experimentation with Principal Component Pursuit (PCP) algorithms.
#'
#' @details
#' The data is simulated as follows:
#'
#' `L <- matrix(runif(n * r), n, r) %*% matrix(runif(r * p), r, p)`
#'
#' `S <- matrix(0, n, p)`
#'
#' `S[sparse_nonzero_idxs] <- 1`
#'
#' `Z <- matrix(rnorm(n * p, sd = sigma), n, p)`
#'
#' `D <- L + S + Z`
#'
#' @param n,p (Optional) A pair of integers specifying the simulated dataset's
#'   number of `n` observations (rows) and `p` variables (columns). By default,
#'   `n = 100`, and `p = 10`.
#' @param r (Optional) An integer specifying the rank of the simulated dataset's
#'   low-rank component. Intuitively, the number of latent patterns governing
#'   the simulated dataset. Must be that `r` `<= min(n, p)`.
#'   By default, `r` = 3.
#' @param sparse_nonzero_idxs (Optional) An integer vector with
#'   `length(sparse_nonzero_idxs) <= n * p` specifying the
#'   indices of the non-zero elements in the sparse component. By default,
#'   `sparse_nonzero_idxs = NULL`, in which case it is defined to be the
#'   vector `seq(1, n * p, n + 1)` (placing sparse noise along the diagonal
#'   of the simulated dataset).
#' @param sigma (Optional) A double specifying the standard deviation of the
#'   dense (Gaussian) noise component `Z`. By default, `sigma = 0.05`.
#' @param seed (Optional) An integer specifying the seed for random number
#'   generation. By default, `seed = 42`.
#'
#' @returns A list containing:
#' * `D`: The observed data matrix, where `D = L + S + Z`.
#' * `L`: The ground truth rank-`r` low-rank matrix.
#' * `S`: The ground truth sparse matrix.
#' * `S`: The ground truth dense (Gaussian) noise matrix.
#'
#' @seealso [sim_na()], [sim_lod()], [impute_matrix()]
#' @examples
#' # rank 3 example
#' data <- sim_data()
#' matrix_rank(data$D)
#' matrix_rank(data$L)
#' # rank 7 example
#' data <- sim_data(n = 1000, p = 25, r = 7)
#' matrix_rank(data$D)
#' matrix_rank(data$L)
#' @export
#' @importFrom stats rnorm runif
sim_data <- function(n = 100, p = 10, r = 3, sparse_nonzero_idxs = NULL, sigma = 0.05, seed = 42) {
  checkmate::qassert(n, rules = "X1[1,)")
  checkmate::qassert(p, rules = "X1[1,)")
  checkmate::qassert(r, rules = "X1[1,)")
  checkmate::assert_true(r <= min(n, p))
  if (is.null(sparse_nonzero_idxs)) sparse_nonzero_idxs <- seq(1, n * p, n + 1)
  checkmate::assert_integerish(sparse_nonzero_idxs)
  checkmate::assert_true(length(sparse_nonzero_idxs) <= n * p)
  checkmate::qassert(sigma, rules = "N1[0,)")
  checkmate::qassert(seed, rules = "X1")
  set.seed(seed)
  L <- matrix(runif(n * r), n, r) %*% matrix(runif(r * p), r, p)
  S <- matrix(0, n, p)
  S[sparse_nonzero_idxs] <- 1
  Z <- matrix(rnorm(n * p, sd = sigma), n, p)
  D <- L + S + Z
  list(D = D, L = L, S = S, Z = Z)
}
