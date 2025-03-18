#' Retrieve default PCP parameter settings for given matrix
#'
#' @description
#' `get_pcp_defaults()` calculates "default" PCP parameter settings `lambda`,
#' `mu` (used in [root_pcp()]), and `eta` (used in [rrmc()]) for a given data
#' matrix `D`.
#'
#' The "default" values of `lambda` and `mu` offer _theoretical_ guarantees
#' of optimal estimation performance.
#' [Candès et al. (2011)](https://doi.org/10.1145/1970392.1970395) obtained
#' the guarantee for `lambda`, while
#' [Zhang et al. (2021)](https://proceedings.neurips.cc/paper/2021/hash/f65854da4622c1f1ad4ffeb361d7703c-Abstract.html)
#' obtained the result for `mu`. _It has not yet been proven whether or
#' not `eta` enjoys similar properties._
#'
#'  _In practice_ it is common to find different optimal parameter values
#' after tuning these parameters in a grid search. Therefore, **it is
#' recommended to use these defaults primarily to help define a reasonable
#' initial parameter search space to pass into [grid_search_cv()].**
#'
#' @section The intuition behind PCP parameters:
#' [root_pcp()]'s objective function is given by:
#' \deqn{\min_{L, S} ||L||_* + \lambda ||S||_1 + \mu ||L + S - D||_F}
#' * `lambda` controls the sparsity of [root_pcp()]'s output `S` matrix;
#'   larger values of `lambda` penalize non-zero entries in `S` more
#'   stringently, driving the recovery of sparser `S` matrices. Therefore,
#'   if you a priori expect few outlying events in your model, you might
#'   expect a grid search to recover relatively larger `lambda` values, and
#'   vice-versa.
#' * `mu` adjusts [root_pcp()]'s sensitivity to noise; larger values of `mu`
#'   penalize errors between the predicted model and the observed data (i.e.
#'   noise), more severely. Environmental data subject to higher noise levels
#'   therefore require a [root_pcp()] model equipped with smaller `mu` values
#'   (since higher noise means a greater discrepancy between the observed
#'   mixture and the true underlying low-rank and sparse model). In virtually
#'   noise-free settings (e.g. simulations), larger values of `mu` would be
#'   appropriate.
#'
#' [rrmc()]'s objective function is given by:
#' \deqn{\min_{L, S} I_{rank(L) \leq r} + \eta ||S||_0 + ||L + S - D||_F^2}
#' * `eta` controls the sparsity of [rrmc()]'s output `S` matrix, just as
#'   `lambda` does for [root_pcp()]. Because there are no other parameters
#'   scaling the noise term, `eta` can be thought of as a ratio between
#'   [root_pcp()]'s `lambda` and `mu`: Larger values of `eta` will place a
#'   greater emphasis on penalizing the non-zero entries in `S` over penalizing
#'   the errors between the predicted and observed data (the dense noise `Z`).
#'
#' @section The calculation of the "default" PCP parameters:
#' * `lambda` is calculated as \eqn{\lambda = 1 / \sqrt{\max(n, p)},} where
#'   \eqn{n} and \eqn{p} are the dimensions of the input matrix
#'   \eqn{D_{n \times p}}
#'   [[Candès et al. (2011)](https://doi.org/10.1145/1970392.1970395)].
#' * `mu` is calculated as \eqn{\mu = \sqrt{\frac{\min(n, p)}{2}},} where
#'   \eqn{n} and \eqn{p} are as above
#'   [[Zhang et al. (2021)](https://proceedings.neurips.cc/paper/2021/hash/f65854da4622c1f1ad4ffeb361d7703c-Abstract.html)].
#' * `eta` is simply \eqn{\eta = \frac{\lambda}{\mu}}.
#'
#' @param D The input data matrix.
#'
#' @returns A list containing:
#' * `lambda`: The theoretically optimal `lambda` value used in [root_pcp()].
#' * `mu`: The theoretically optimal `mu` value used in [root_pcp()].
#' * `eta`: The default `eta` value used in [rrmc()].
#'
#' @seealso [grid_search_cv()]
#' @examples
#' # Examine the queens PM2.5 data
#' queens
#' # Get rid of the Date column
#' D <- queens[, 2:ncol(queens)]
#' # Get default PCP parameters
#' default_params <- get_pcp_defaults(D)
#' # Use default parameters to define parameter search space
#' scaling_factors <- sort(c(10^seq(-2, 4, 1), 2 * 10^seq(-2, 4, 1)))
#' etas_to_grid_search <- default_params$eta * scaling_factors
#' etas_to_grid_search
#' @references Candès, Emmanuel J., Xiaodong Li, Yi Ma, and John Wright.
#'   "Robust principal component analysis?." Journal of the ACM (JACM)
#'   58, no. 3 (2011): 1-37. [available
#'   [here](https://doi.org/10.1145/1970392.1970395)]
#' @references Zhang, Junhui, Jingkai Yan, and John Wright.
#'   "Square root principal component pursuit: tuning-free noisy robust matrix
#'   recovery." Advances in Neural Information Processing Systems 34 (2021):
#'   29464-29475. [available
#'   [here](https://proceedings.neurips.cc/paper/2021/hash/f65854da4622c1f1ad4ffeb361d7703c-Abstract.html)]
#' @export
get_pcp_defaults <- function(D) {
  n <- nrow(D)
  p <- ncol(D)
  lambda <- 1 / sqrt(max(n, p))
  mu <- sqrt(min(n, p) / 2)
  eta <- lambda / mu
  list(lambda = lambda, mu = mu, eta = eta)
}

#' Simulate data
#'
#' `sim_data()` generates a simulated dataset `D = L + S + Z` for
#' experimentation with Principal Component Pursuit (PCP) algorithms.
#'
#' @param n,p (Optional) A pair of integers specifying the simulated dataset's
#'   number of `n` observations (rows) and `p` variables (columns). By default,
#'   `n = 100`, and `p = 10`.
#' @param r (Optional) An integer specifying the rank of the simulated dataset's
#'   low-rank component. Intuitively, the number of latent patterns governing
#'   the simulated dataset. By default, `r` = 3.
#' @param sparse_nonzero_idxs (Optional) An integer vector specifying the
#'   indices of the non-zero elements in the sparse component. By default,
#'   `sparse_nonzero_idxs = NULL`, in which case it is defined to be the
#'   vector `c(1, n * p)` (placing sparse noise in the very first and last
#'   entries of the simulated dataset).
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
#' @seealso [corrupt_mat_randomly()]
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
  set.seed(seed)
  if (is.null(sparse_nonzero_idxs)) sparse_nonzero_idxs <- c(1, n * p)
  L <- matrix(runif(n * r), n, r) %*% matrix(runif(r * p), r, p)
  S <- matrix(0, n, p)
  S[sparse_nonzero_idxs] <- 1
  Z <- matrix(rnorm(n * p, sd = sigma), n, p)
  D <- L + S + Z
  list(D = D, L = L, S = S, Z = Z)
}

#' Compute singular values of given matrix
#'
#' @description
#' `sing()` calculates the singular values of a given data matrix `D`. This is
#' done with a call to [svd()], and is included in `pcpr` to enable the quick
#' characterization of a data matrix's raw low-rank structure, to help decide
#' whether [rrmc()] or [root_pcp()] is the more appropriate PCP algorithm to
#' employ in conjunction with `D`.
#'
#' Experimentally, the [rrmc()] approach to PCP has best been able to handle
#' those datasets that are governed by complex underlying patterns characterized
#' by slowly decaying singular values, such as EH data. For observed data with a
#' well-defined low rank structure (rapidly decaying singular values),
#' [root_pcp()] may offer a better model estimate.
#'
#' @param D The input data matrix.
#'
#' @returns A numeric vector containing the singular values of `D`.
#' @examples
#' data <- sim_data()
#' sing(data$D)
#' # could plot the singular values for visual inspection with e.g.
#' # plot(sing(data$D), type = 'b')
#' @references "Singular value decomposition" [Wikipedia article](https://en.wikipedia.org/wiki/Singular_value_decomposition).
#' @export
sing <- function(D) {
  svd(D)$d
}