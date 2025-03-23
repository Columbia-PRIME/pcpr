#' Retrieve default PCP parameter settings for given matrix
#'
#' @description
#' `get_pcp_defaults()` calculates "default" PCP parameter settings `lambda`,
#' `mu` (used in [root_pcp()]), and `eta` (used in [rrmc()]) for a given data
#' matrix `D`.
#'
#' The "default" values of `lambda` and `mu` offer _theoretical_ guarantees
#' of optimal estimation performance. Candès et al. (2011) obtained the
#' guarantee for `lambda`, while
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
#'   \eqn{D_{n \times p}} Candès et al. (2011).
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
#' D <- as.matrix(queens[, 2:ncol(queens)])
#' # Get default PCP parameters
#' default_params <- get_pcp_defaults(D)
#' # Use default parameters to define parameter search space
#' scaling_factors <- sort(c(10^seq(-2, 4, 1), 2 * 10^seq(-2, 4, 1)))
#' etas_to_grid_search <- default_params$eta * scaling_factors
#' etas_to_grid_search
#' @references Candès, Emmanuel J., Xiaodong Li, Yi Ma, and John Wright.
#'   "Robust principal component analysis?." Journal of the ACM (JACM)
#'   58, no. 3 (2011): 1-37.
#' @references Zhang, Junhui, Jingkai Yan, and John Wright.
#'   "Square root principal component pursuit: tuning-free noisy robust matrix
#'   recovery." Advances in Neural Information Processing Systems 34 (2021):
#'   29464-29475. [available
#'   [here](https://proceedings.neurips.cc/paper/2021/hash/f65854da4622c1f1ad4ffeb361d7703c-Abstract.html)]
#' @export
get_pcp_defaults <- function(D) {
  checkmate::assert_matrix(D)
  n <- nrow(D)
  p <- ncol(D)
  lambda <- 1 / sqrt(max(n, p))
  mu <- sqrt(min(n, p) / 2)
  eta <- lambda / mu
  list(lambda = lambda, mu = mu, eta = eta)
}
