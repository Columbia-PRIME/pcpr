#' Square root principal component pursuit (convex PCP)
#'
#' @description
#' `root_pcp()` implements the convex PCP algorithm "Square root principal
#' component pursuit" as described in
#' [Zhang et al. (2021)](https://proceedings.neurips.cc/paper/2021/hash/f65854da4622c1f1ad4ffeb361d7703c-Abstract.html)
#' , outfitted with environmental health (EH)-specific extensions as described
#' in Gibson et al. (2022).
#'
#' Given an observed data matrix `D`, and regularization parameters `lambda` and
#' `mu`, `root_pcp()` aims to find the best low-rank and sparse estimates `L`
#' and `S`. The `L` matrix encodes latent patterns that govern the observed
#' data. The `S` matrix captures any extreme events in the data unexplained by
#' the underlying patterns in `L`.
#'
#' Being convex, `root_pcp()` determines the rank `r`, or number of latent
#' patterns in the data, autonomously during it's optimization. As such, the
#' user does not need to specify the desired rank `r` of the output `L` matrix
#' as in the non-convex PCP model [rrmc()].
#'
#' Experimentally, the `root_pcp()` approach to PCP modeling has best been able
#' to handle those datasets that are governed by well-defined underlying
#' patterns, characterized by quickly decaying singular values. This is typical
#' of imaging and video data, but uncommon for EH data. For observed data with a
#' complex low rank structure (slowly decaying singular values), like EH data,
#' [rrmc()] may offer a better model estimate.
#'
#' Three EH-specific extensions are currently supported by `root_pcp()`:
#' 1. The model can handle missing values in the input data matrix `D`;
#' 2. The model can also handle measurements that fall below the limit of
#'    detection (LOD), if provided `LOD` information by the user; and
#' 3. The model is also equipped with an optional non-negativity constraint
#'    on the low-rank `L` matrix, ensuring that all output values in `L` are
#'    \eqn{> 0}.
#'
#' @section The objective function:
#' `root_pcp()` optimizes the following objective function:
#' \deqn{\min_{L, S} ||L||_* + \lambda ||S||_1 + \mu ||L + S - D||_F}
#' The first term is the nuclear norm of the `L` matrix, incentivizing `L` to be
#' low-rank. The second term is the \eqn{\ell_1} norm of the `S` matrix,
#' encouraging `S` to be sparse. The third term is the Frobenius norm
#' applied to the model's noise, ensuring that the estimated low-rank and sparse
#' models `L` and `S` together have high fidelity to the observed data `D`.
#' The objective is not smooth nor differentiable, however it is convex and
#' separable. As such, it is optimized using the Alternating Direction
#' Method of Multipliers (ADMM) algorithm Boyd et al. (2011), Gao et al. (2020).
#'
#' @section The `lambda` and `mu` parameters:
#' * `lambda` controls the sparsity of `root_pcp()`'s output `S` matrix;
#'   larger values of `lambda` penalize non-zero entries in `S` more
#'   stringently, driving the recovery of sparser `S` matrices. Therefore,
#'   if you a priori expect few outlying events in your model, you might
#'   expect a grid search to recover relatively larger `lambda` values, and
#'   vice-versa.
#' * `mu` adjusts `root_pcp()`'s sensitivity to noise; larger values of `mu`
#'   penalize errors between the predicted model and the observed data (i.e.
#'   noise), more severely. Environmental data subject to higher noise levels
#'   therefore require a `root_pcp()` model equipped with smaller `mu` values
#'   (since higher noise means a greater discrepancy between the observed
#'   mixture and the true underlying low-rank and sparse model). In virtually
#'   noise-free settings (e.g. simulations), larger values of `mu` would be
#'   appropriate.
#'
#' The default values of `lambda` and `mu` offer _theoretical_ guarantees
#' of optimal estimation performance, and stable recovery of `L` and `S`. By
#' "stable", we mean `root_pcp()`'s reconstruction error is, in the worst case,
#' proportional to the magnitude of the noise corrupting the observed data
#' (\eqn{||Z||_F}), often outperforming this upper bound.
#' Candès et al. (2011) obtained the guarantee for `lambda`, while
#' [Zhang et al. (2021)](https://proceedings.neurips.cc/paper/2021/hash/f65854da4622c1f1ad4ffeb361d7703c-Abstract.html)
#' obtained the result for `mu`.
#'
#' @section Environmental health specific extensions:
#' We refer interested readers to
#' Gibson et al. (2022) for the complete details regarding the EH-specific
#' extensions.
#'
#' **Missing value functionality:** PCP assumes that the same data generating
#' mechanisms govern both the missing and the observed entries in `D`. Because
#' PCP primarily seeks accurate estimation of _patterns_ rather than
#' individual _observations_, this assumption is reasonable, but in some edge
#' cases may not always be justified. Missing values in `D` are therefore
#' reconstructed in the recovered low-rank `L` matrix according to the
#' underlying patterns in `L`. There are three corollaries to keep in mind
#' regarding the quality of recovered missing observations:
#' 1. Recovery of missing entries in `D` relies on accurate estimation of
#'    `L`;
#' 2. The fewer observations there are in `D`, the harder it is to accurately
#'    reconstruct `L` (therefore estimation of _both_ unobserved _and_ observed
#'    measurements in `L` degrades); and
#' 3. Greater proportions of missingness in `D` artifically drive up the
#'    sparsity of the estimated `S` matrix. This is because it is not possible
#'    to recover a sparse event in `S` when the corresponding entry in `D` is
#'    unobserved. By definition, sparse events in `S` cannot be explained by
#'    the consistent patterns in `L`. Practically, if 20% of the entries in `D`
#'    are missing, then at least 20% of the entries in `S` will be 0.
#'
#' **Handling measurements below the limit of detection:** When equipped with
#' LOD information, PCP treats any estimations of values known to be below the
#' LOD as equally valid if their approximations fall between 0 and the LOD. Over
#' the course of optimization, observations below the LOD are pushed into this
#' known range \eqn{[0, LOD]} using penalties from above and below: should a
#' \eqn{< LOD} estimate be \eqn{< 0}, it is stringently penalized, since
#' measured observations cannot be negative. On the other hand, if a \eqn{< LOD}
#' estimate is \eqn{>} the LOD, it is also heavily penalized: less so than when
#' \eqn{< 0}, but more so than observations known to be above the LOD, because
#' we have prior information that these observations must be below LOD.
#' Observations known to be above the LOD are penalized as usual, using the
#' Frobenius norm in the above objective function.
#'
#' Gibson et al. (2022) demonstrates that
#' in experimental settings with up to 50% of the data corrupted below the LOD,
#' PCP with the LOD extension boasts superior accuracy of recovered `L` models
#' compared to PCA coupled with \eqn{LOD / \sqrt{2}} imputation. PCP even
#' outperforms PCA in low-noise scenarios with as much as 75% of the data
#' corrupted below the LOD. The few situations in which PCA bettered PCP were
#' those pathological cases in which `D` was characterized by extreme noise and
#' huge proportions (i.e., 75%) of observations falling below the LOD.
#'
#' **The non-negativity constraint on `L`:** To enhance interpretability of
#' PCP-rendered solutions, there is an optional non-negativity constraint
#' that can be imposed on the `L` matrix to ensure all estimated values
#' within it are \eqn{\geq 0}. This prevents researchers from having to deal
#' with negative observation values and questions surrounding their meaning
#' and utility. Non-negative `L` models also allow for seamless use of methods
#' such as non-negative matrix factorization to extract non-negative patterns.
#' The non-negativity constraint is incorporated in the ADMM splitting technique
#' via the introduction of an additional optimization variable and corresponding
#' constraint.
#' @inheritParams rrmc
#' @param lambda,mu (Optional) A pair of doubles each in the range `[0, Inf)`
#'   regularizing `S` and `L`. `lambda` controls the sparsity of the output
#'   `S` matrix; larger values penalize non-zero entries in `S` more
#'   stringently, driving the recovery of sparser `S` matrices. `mu` adjusts the
#'    model's sensitivity to noise; larger values will penalize errors between
#'   the predicted model and the observed data more severely. It is highly
#'   recommended the user tunes both of these parameters using
#'   [grid_search_cv()] for each unique data matrix `D`. By default, both
#'   `lambda` and `mu` are `NULL`, in which case the theoretically optimal
#'   values are used, calculated according to [get_pcp_defaults()].
#' @param non_negative (Optional) A logical indicating whether or not the
#'   non-negativity constraint should be used to constrain the output `L`
#'   matrix to have all entries \eqn{\geq 0}. By default, `non_negative = TRUE`.
#' @param max_iter (Optional) An integer specifying the maximum number of
#'   iterations to allow PCP before giving up on meeting PCP's convergence
#'   criteria. By default, `max_iter = 10000`, suitable for most problems.
#' @param verbose (Optional) A logical indicating whether or not to print
#'   information in real time over the course of PCP's optimization. By
#'   default, `verbose = FALSE`.
#'
#' @returns A list containing:
#'   * `L`: The rank-`r` low-rank matrix encoding the `r`-many latent patterns
#'     governing the observed input data matrix `D`. `dim(L)` will be the same
#'     as `dim(D)`. To explicitly obtain the underlying patterns, `L` can be
#'     used as the input to any matrix factorization technique of choice, e.g.
#'     PCA, factor analysis, or non-negative matrix factorization.
#'   * `S`: The sparse matrix containing the rare outlying or extreme
#'     observations in `D` that are not explained by the underlying patterns in
#'     the corresponding `L` matrix. `dim(S)` will be the same as `dim(D)`.
#'     Most entries in `S` are `0`, while non-zero entries identify the extreme
#'     outlying observations in `D`.
#'   * `num_iter`: The number of iterations taken to reach convergence. If
#'     `num_iter == max_iter` then `root_pcp()` did not converge.
#'   * `objective`: A vector containing the values of `root_pcp()`'s objective
#'     function over the course of optimization.
#'   * `converged`: A boolean indicating whether the convergence criteria were
#'     met before `max_iter` was reached.
#'
#' @seealso [rrmc()]
#' @examples
#' #### -------Simple simulated PCP problem-------####
#' # First we will simulate a simple dataset with the sim_data() function.
#' # The dataset will be a 100x10 matrix comprised of:
#' # 1. A rank-2 component as the ground truth L matrix;
#' # 2. A ground truth sparse component S w/outliers along the diagonal; and
#' # 3. A dense Gaussian noise component
#' data <- sim_data(r = 2, sigma = 0.1)
#' # Normally we would conduct grid search to tune lambda and mu. But, to keep
#' # the example short, we will just use best parameters found in the below grid
#' # search example:
#' \dontrun{
#' lambda_0 <- get_pcp_defaults(data$D)$lambda
#' mu_0 <- get_pcp_defaults(data$D)$mu
#' lambdas <- lambda_0 + seq(-0.05, 0.2, 0.025)
#' mus <- mu_0 + seq(-1, 1, 0.3)
#' params <- expand.grid(lambdas, mus)
#' names(params) <- c("lambda", "mu")
#' gs <- grid_search_cv(data$D, root_pcp, params)
#' dplyr::arrange(gs$summary_stats, rel_err)
#' }
#' # The gs found the best parameters to be lambda = 0.225 and mu = 3.04
#' pcp_model <- root_pcp(data$D, lambda = 0.225, mu = 3.04)
#' data.frame(
#'   "Estimated_L_rank" = matrix_rank(pcp_model$L, 5e-2),
#'   "Observed_relative_error" = norm(data$L - data$D, "F") / norm(data$L, "F"),
#'   "PCA_error" = norm(data$L - proj_rank_r(data$D, r = 2), "F") / norm(data$L, "F"),
#'   "PCP_L_error" = norm(data$L - pcp_model$L, "F") / norm(data$L, "F"),
#'   "PCP_S_error" = norm(data$S - pcp_model$S, "F") / norm(data$S, "F")
#' )
#' # Results:
#' # PCP found a rank 2 solution!
#' # PCP outperformed PCA in it's recovery of the L matrix (even though we let
#' # PCA "cheat" by telling PCA it was looking for a rank 2 solution)!
#' # PCP successfully isolated the outlying events in S!
#' @references Zhang, Junhui, Jingkai Yan, and John Wright.
#'   "Square root principal component pursuit: tuning-free noisy robust matrix
#'   recovery." Advances in Neural Information Processing Systems 34 (2021):
#'   29464-29475. [available
#'   [here](https://proceedings.neurips.cc/paper/2021/hash/f65854da4622c1f1ad4ffeb361d7703c-Abstract.html)]
#' @references Gibson, Elizabeth A., Junhui Zhang, Jingkai Yan, Lawrence
#'   Chillrud, Jaime Benavides, Yanelli Nunez, Julie B. Herbstman, Jeff
#'   Goldsmith, John Wright, and Marianthi-Anna Kioumourtzoglou.
#'   "Principal component pursuit for pattern identification in
#'   environmental mixtures." Environmental Health Perspectives 130, no.
#'   11 (2022): 117008.
#' @references Boyd, Stephen, Neal Parikh, Eric Chu, Borja Peleato, and Jonathan
#'   Eckstein. "Distributed optimization and statistical learning via the
#'   alternating direction method of multipliers." Foundations and Trends in
#'   Machine learning 3, no. 1 (2011): 1-122.
#' @references Gao, Wenbo, Donald Goldfarb, and Frank E. Curtis. "ADMM for
#'   multiaffine constrained optimization." Optimization Methods and Software
#'   35, no. 2 (2020): 257-303.
#' @references Candès, Emmanuel J., Xiaodong Li, Yi Ma, and John Wright.
#'   "Robust principal component analysis?." Journal of the ACM (JACM)
#'   58, no. 3 (2011): 1-37.
#' @export
root_pcp <- function(D, lambda = NULL, mu = NULL, LOD = -Inf, non_negative = TRUE, max_iter = 10000, verbose = FALSE) {
  # 1. Initialize variables, argument assertions:
  # 1a. Matrix dimensions:
  checkmate::assert_matrix(D)
  n <- nrow(D)
  p <- ncol(D)
  if (is.null(lambda)) lambda <- 1 / sqrt(max(n, p))
  if (is.null(mu)) mu <- sqrt(min(n, p) / 2)
  checkmate::qassert(lambda, rules = "N1[0,)")
  checkmate::qassert(mu, rules = "N1[0,)")
  checkmate::qassert(LOD, rules = "n+")
  if (checkmate::test_numeric(LOD, min.len = 2) && !checkmate::test_matrix(LOD)) {
    checkmate::assert_true(ncol(D) == length(LOD))
  }
  if (checkmate::test_matrix(LOD)) {
    checkmate::assert_true(all(dim(D) == dim(LOD)))
  }
  checkmate::qassert(non_negative, rules = "B1")
  checkmate::qassert(max_iter, rules = "X1[1,)")
  checkmate::qassert(verbose, rules = "B1")
  # 1b. Process LOD input (convert LOD vector to matrix if needed, so it multiplies correctly):
  if (is.vector(LOD) && length(LOD) != 1) LOD <- matrix(LOD, nrow = n, ncol = p, byrow = TRUE)
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
  L1 <- matrix(0, nrow = n, ncol = p)
  L2 <- matrix(0, nrow = n, ncol = p)
  L3 <- matrix(0, nrow = n, ncol = p)
  # 1g. S matrix primal variables:
  S1 <- matrix(0, nrow = n, ncol = p)
  S2 <- matrix(0, nrow = n, ncol = p)
  # 1h. Noise/error matrix primal variable:
  Z <- matrix(0, nrow = n, ncol = p)
  # 1i. Dual variables (noise matrix; Y4 carries non-negativity constraint in dual):
  Y1 <- matrix(0, nrow = n, ncol = p)
  Y2 <- matrix(0, nrow = n, ncol = p)
  Y3 <- matrix(0, nrow = n, ncol = p)
  Y4 <- matrix(0, nrow = n, ncol = p)
  # 1j. Flag identifying convergence; vector of objective values for user
  converged <- FALSE
  objective <- vector("numeric", length = max_iter)
  # 2. ADMM-splitting iterations:
  for (i in 1:max_iter) {
    # 2a. Update 1st primal variables (L1, S1):
    nuc <- prox_nuclear((L2 + L3 - (Y1 + Y4) / rho) / 2, c = 1 / 2 / rho)
    L1 <- nuc[[1]]
    L1_nuclear_norm <- nuc[[2]]
    S1 <- prox_l1(S2 - Y2 / rho, c = lambda / rho)
    # 2b. Update 2nd primal variables (L2, L3, S2):
    L2_old <- L2
    S2_old <- S2
    L3_old <- L3
    temp <- L2 + S2 - Y3 / rho
    temp_D <- D * mask_above_LOD + temp * (mask_below_LOD & (temp >= 0) & (temp <= LOD)) + ifelse(mask_below_LOD & (temp > LOD), yes = LOD, no = 0)
    Z <- prox_frobenius(temp - temp_D, c = mu / rho) + temp_D
    if (non_negative) L3 <- pmax(L1 + Y4 / rho, 0)
    L2_obs <- Omega / 3 * (2 * L1 - S1 + Z + (2 * Y1 - Y2 + Y3) / rho)
    L2_unobs <- (1 - Omega) * (L1 + Y1 / rho)
    L2 <- L2_obs + L2_unobs
    S2_obs <- Omega / 3 * (2 * S1 - L1 + Z + (2 * Y2 - Y1 + Y3) / rho)
    S2_unobs <- (1 - Omega) * (S1 + Y2 / rho)
    S2 <- S2_obs + S2_unobs
    # 2c. Update dual variables (Y1, Y2, Y3, Y4):
    Y1 <- Y1 + rho * (L1 - L2)
    Y2 <- Y2 + rho * (S1 - S2)
    Y3 <- Y3 + rho * Omega * (Z - (L2 + S2))
    if (non_negative) Y4 <- Y4 + rho * (L1 - L3)
    # 2d. Calculate primal & dual residuals:
    res_primal <- sqrt(norm(L1 - L2, type = "F")^2 + non_negative * norm(L1 - L3, type = "F")^2 + norm(S1 - S2, type = "F")^2 + norm(Omega * (Z - L2 - S2), type = "F")^2)
    res_dual <- rho * sqrt(norm(L2 + L3 - L2_old - L3_old, type = "F")^2 + norm(S2 - S2_old, type = "F")^2 + norm(Omega * (L2 - L2_old + S2 - S2_old), type = "F")^2)
    # 2e. Update rho:
    if (res_primal > 10 * res_dual) {
      rho <- 2 * rho
    } else if (res_dual > 10 * res_primal) {
      rho <- rho / 2
    }
    # 2f. Calculate objective:
    objective[i] <- L1_nuclear_norm + lambda * sum(abs(S1)) + mu * sqrt(loss_lod(D, L2 + S2, LOD))
    if (verbose) cat(paste0("Iteration: ", i, " Obj: ", round(objective[i], 5)))
    # 2g. Update & check stopping criteria:
    thresh_primal <- eps_abs * sqrt((5 + non_negative) * n * p) + eps_rel * max(sqrt((1 + non_negative) * norm(L1, type = "F")^2 + norm(S1, type = "F")^2 + norm(Z, type = "F")^2), sqrt(norm(L2, type = "F")^2 + norm(L3, type = "F")^2 + norm(S2, type = "F")^2 + norm(Omega * (L2 + S2), type = "F")^2))
    thresh_dual <- eps_abs * sqrt((3 + non_negative) * n * p) + eps_rel * sqrt(norm(Y1 + Y4, type = "F")^2 + norm(Y2, type = "F")^2 + norm(Y3, type = "F")^2)
    if (res_primal < thresh_primal && res_dual < thresh_dual) {
      converged <- TRUE
      if (verbose) cat(paste0("Converged in ", i, " iterations."))
      break
    }
  }
  # 3. Wrap up & return:
  if (!converged && verbose) warning(paste0("Maximum iterations reached (", max_iter, "). PCP did not converge."))
  L_final <- (L1 + L2 + L3) / (2 + non_negative)
  S_final <- (S1 + S2) / 2
  if (non_negative) L_final[L_final < 0] <- 0
  list(L = L_final, S = S_final, num_iter = i, objective = objective, converged = converged)
}

#' Loss function for PCP's limit of detection penalty
#'
#' `loss_lod()` includes the LOD-specific penalty terms to compute the loss on
#' the squared error term `L + S - D`.
#'
#' @param D The original data matrix.
#' @param X The predicted value of `L + S` at the current iteration.
#' @param LOD The LOD as a matrix, where `dim(LOD) == dim(D)`.
#'
#' @returns Scalar value used to calculate loss in objective function.
#'
#' @seealso [root_pcp()]
#' @keywords internal
#' @noRd
loss_lod <- function(D, X, LOD) {
  X_lod <- (X - D) * (D >= 0) +
    (X - LOD) * ((D < 0) & (X > LOD)) +
    X * ((D < 0) & (X < 0))
  sum(X_lod^2) / 2
}
