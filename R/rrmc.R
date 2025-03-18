#' Rank-based robust matrix completion (non-convex PCP)
#'
#' @description
#' `rrmc()` implements the non-convex PCP algorithm "Rank-based robust matrix
#' completion" as described in
#' [Cherapanamjeri et al. (2017)](https://proceedings.mlr.press/v70/cherapanamjeri17a.html)
#' (see Algorithm 3), outfitted with environmental health (EH)-specific
#' extensions as described in
#' [Gibson et al. (2022)](https://doi.org/10.1289/EHP10479).
#'
#' Given an observed data matrix `D`, maximum rank to search up to `r`, and
#' regularization parameter `eta`, `rrmc()` seeks to find the best low-rank
#' and sparse estimates `L` and `S` using an incremental rank-based strategy.
#' The `L` matrix encodes latent patterns that govern the observed data.
#' The `S` matrix captures any extreme events in the data unexplained by the
#' underlying patterns in `L`.
#'
#' `rrmc()`'s incremental rank-based strategy first estimates a rank-`1` model
#' \eqn{(L^(1), S^(1))}, before using the rank-`1` model as the initialization
#' point to then construct a rank-`2` model \eqn{(L^(2), S^(2))}, and so on,
#' until the desired rank-`r` model \eqn{(L^(r), S^(r))} is recovered. All
#' models from ranks `1` through `r` are returned by `rrmc()` in this way.
#'
#' Experimentally, the `rrmc()` approach to PCP has best been able to handle
#' those datasets that are governed by complex underlying patterns characterized
#' by slowly decaying singular values, such as EH data. For observed data with a
#' well-defined low rank structure (rapidly decaying singular values),
#' [root_pcp()] may offer a better model estimate.
#'
#' Two EH-specific extensions are currently supported by `rrmc()`:
#' 1. The model can handle missing values in the input data matrix `D`; and
#' 2. The model can also handle measurements that fall below the limit of
#'    detection (LOD), if provided `LOD` information by the user.
#'
#' Support for a non-negativity constraint on `rrmc()`'s output will be added in
#' a future release of `pcpr`.
#'
#' @section The objective function:
#' `rrmc()` implicitly optimizes the following objective function:
#' \deqn{\min_{L, S} I_{rank(L) \leq r} + \eta ||S||_0 + ||L + S - D||_F^2}
#' The first term is the indicator function checking that the `L` matrix is
#' strictly rank `r` or less, implemented using a rank `r` projection operator
#' [proj_rank_r()]. The second term is the \eqn{\ell_0} norm applied to the `S`
#' matrix to encourage sparsity, and is implemented with the help of an adaptive
#' hard-thresholding operator [hard_threshold()]. The third term is the squared
#' Frobenius norm applied to the model's noise.
#'
#' @section The `eta` parameter:
#' The `eta` parameter scales the sparse penalty applied to `rrmc()`'s output
#' sparse `S` matrix. Larger values of `eta` penalize non-zero entries in `S`
#' more stringently, driving the recovery of sparser `S` matrices.
#'
#' Because there are no other parameters scaling the other terms in `rrmc()`'s
#' objective function, `eta` can intuitively be thought of as the dial that
#' balances the model's sensitivity to extreme events (placed in `S`) and
#' its sensitivity to noise `Z` (captured by the last term in the objective,
#' which measures the discrepancy between the between the predicted model
#' and the observed data). Larger values of `eta` will place a
#' greater emphasis on penalizing the non-zero entries in `S` over penalizing
#' the errors between the predicted and observed data `Z = L + S - D`.
#'
#' @section Environmental health specific extensions:
#' We refer interested readers to
#' [Gibson et al. (2022)](https://doi.org/10.1289/EHP10479) for the complete
#' details regarding the EH-specific extensions.
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
#' [Gibson et al. (2022)](https://doi.org/10.1289/EHP10479) demonstrates that
#' in experimental settings with up to 50% of the data corrupted below the LOD,
#' PCP with the LOD extension boasts superior accuracy of recovered `L` models
#' compared to PCA coupled with \eqn{LOD / \sqrt{2}} imputation. PCP even
#' outperforms PCA in low-noise scenarios with as much as 75% of the data
#' corrupted below the LOD. The few situations in which PCA bettered PCP were
#' those pathological cases in which `D` was characterized by extreme noise and
#' huge proportions (i.e., 75%) of observations falling below the LOD.
#'
#' @param D The input data matrix (can contain `NA` values). Note that PCP will
#'   converge much more quickly when `D` has been standardized in some way (e.g.
#'   scaling columns by their standard deviations, or column-wise min-max
#'   normalization).
#' @param r An integer specifying the maximum rank PCP model to return. All
#'   models from rank `1` through `r` will be returned.
#' @param eta (Optional) A double defining the ratio between the model's
#'   sensitivity to sparse and dense noise. Larger values of `eta` will place a
#'   greater emphasis on penalizing the non-zero entries in `S` over penalizing
#'   dense noise `Z`, i.e. errors between the predicted and observed data
#'   `Z = L + S - D`. It is recommended to tune `eta` using [grid_search_cv()]
#'   for each unique data matrix `D`. By default, `eta = NULL`, in which case
#'   `eta` is retrieved using [get_pcp_defaults()].
#' @param LOD (Optional) The limit of detection (LOD) data. Entries in `D` that
#'   satisfy `D >= LOD` are understood to be above the LOD, otherwise those
#'   entries are treated as below the LOD. `LOD` can be either:
#'   * A double, implying a universal LOD common across all measurements in `D`;
#'   * A vector of length `ncol(D)`, signifying a column-specific LOD, where
#'     each entry in the `LOD` vector corresponds to the LOD for each column in
#'     `D`; or
#'   * A matrix of dimension `dim(D)`, indicating an observation-specific LOD,
#'     where each entry in the `LOD` matrix corresponds to the LOD for each
#'     entry in `D`.
#'
#'   By default, `LOD = -Inf`, indicating there are no known LODs for PCP to
#'   leverage.
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
#'   * `L_list`: A list of the `r`-many `L` matrices recovered over the course
#'     of `rrmc()`'s iterative optimization procedure. The first element in
#'     `L_list` corresponds to the rank-`1` `L` matrix, the second to the
#'     rank-`2` `L` matrix, and so on.
#'   * `S_list`: A list of the `r`-many corresponding `S` matrices recovered
#'     over the course of `rrmc()`'s iterative optimization procedure. The first
#'     element in `S_list` corresponds to the rank-`1` solution's `S` matrix,
#'     the second to the rank-`2` solution's `S` matrix, and so on.
#'   * `objective`: A vector containing the values of `rrmc()`'s objective
#'     function over the course of optimization.
#'
#' @seealso [root_pcp()]
#' @examples
#' ####-------Simple simulated PCP problem-------####
#' # First we will simulate a simple dataset with the sim_data() function.
#' # The dataset will be a 100x10 matrix comprised of:
#' # 1. A rank-3 component as the ground truth L matrix;
#' # 2. A ground truth sparse component S w/outliers in 1st & last entries; and
#' # 3. A dense Gaussian noise component
#' data <- sim_data()
#' # Normally we would conduct grid search to tune eta. But, to keep the example
#' # short, we will just use best parameters from the below grid search example:
#' \dontrun{
#' eta_0 <- get_pcp_defaults(data$D)$eta
#' eta_grid <- data.frame("eta" = sort(c(0.1 * eta_0, eta_0 * seq(1, 10, 2))), "r" = 7)
#' gs <- grid_search_cv(data$D, rrmc, eta_grid)
#' dplyr::arrange(gs$summary_stats, rel_err)
#' }
#' # The gs found the best rank to be 3, and the best eta to be 0.3 or 0.4, so
#' # we will split the difference and use an eta of 0.35
#' pcp_model <- rrmc(data$D, r = 3, eta = 0.35)
#' data.frame(
#'   "Observed_relative_error" = norm(data$L - data$D, "F") / norm(data$L, "F"),
#'   "PCA_error" = norm(data$L - proj_rank_r(data$D, r = 3), "F") / norm(data$L, "F"),
#'   "PCP_L_error" = norm(data$L - pcp_model$L, "F") / norm(data$L, "F"),
#'   "PCP_S_error" = norm(data$S - pcp_model$S, "F") / norm(data$S, "F")
#' )
#' # Results:
#' # The grid search correctly found the rank (3) of the ground truth L matrix!
#' # PCP outperformed PCA in it's recovery of the L matrix (even though we let
#' # PCA "cheat" by telling PCA it was looking for a rank 3 solution)!
#' # PCP successfully isolated the outlying event in S!
#' @references Cherapanamjeri, Yeshwanth, Kartik Gupta, and Prateek Jain.
#'   "Nearly optimal robust matrix completion."
#'   International Conference on Machine Learning. PMLR, 2017. [available
#'   [here](https://proceedings.mlr.press/v70/cherapanamjeri17a.html)]
#' @references Gibson, Elizabeth A., Junhui Zhang, Jingkai Yan, Lawrence
#'   Chillrud, Jaime Benavides, Yanelli Nunez, Julie B. Herbstman, Jeff
#'   Goldsmith, John Wright, and Marianthi-Anna Kioumourtzoglou.
#'   "Principal component pursuit for pattern identification in
#'   environmental mixtures." Environmental Health Perspectives 130, no.
#'   11 (2022): 117008. [available [here](https://doi.org/10.1289/EHP10479)]
#' @export
rrmc <- function(D, r, eta = NULL, LOD = -Inf) {
  # 1. Initialize variables:
  # 1a. Matrix dimensions:
  n <- nrow(D)
  p <- ncol(D)
  if (is.null(eta)) eta <- (1 / sqrt(max(n, p))) / (sqrt(min(n, p) / 2))
  # 1b. Process LOD input (convert LOD vector to matrix if needed, so it multiplies correctly):
  if (any(class(LOD) == "list")) LOD <- unlist(LOD)
  if (is.vector(LOD) && length(LOD) != 1) LOD <- matrix(LOD, n, p, TRUE)
  # 1c. Omega = mask of observed values, for handling missingness later:
  Omega <- !is.na(D)
  D[!Omega] <- 0
  # 1d. mask_above/below_LOD = entries that are observed and above/below LOD:
  mask_above_LOD <- Omega & (D >= LOD)
  mask_below_LOD <- Omega & (D < LOD)
  # 1e. Hard-coded optimization parameters:
  epsilon <- 1e-05 * norm(D, "F")
  temp <- svd(D)$d[1]
  t <- 10 * log(40 * r * n * temp / epsilon)
  zeta <- eta * temp
  # 1f. Instantiate variables to return to user:
  L <- matrix(0, nrow = n, ncol = p)
  L_list <- list()
  S_list <- list()
  objective <- vector("numeric", r * t)
  objective_index <- 1
  # 2. Iterative optimization procedure:
  for (k in 1:r) {
    for (i in 0:t) {
      S <- hard_threshold(Omega * (D - L), zeta)
      D_hat <- L + S # approximation to the data
      # calculate the gradient of the LOD penalty:
      grad_LOD_penalty <- mask_above_LOD * (D_hat - D) + # observed and above LOD
        mask_below_LOD * (D_hat < 0) * D_hat + # observed, below LOD, less than 0
        ifelse(mask_below_LOD, (D_hat > LOD) * (D_hat - LOD), 0) # observed, below LOD, current appx bigger than LOD
      # gradient step:
      M_i <- L - (n * p / sum(Omega)) * grad_LOD_penalty
      L <- proj_rank_r(M_i, k)
      temp <- svd(M_i)$d
      zeta <- eta * (temp[k + 1] + 0.5^(i - 2) * temp[k])
      # calculate objective function value:
      objective[objective_index] <- norm(D - L - S, "F")^2 + eta * sum(S != 0)
      objective_index <- objective_index + 1
    }
    L_list[[k]] <- L
    S_list[[k]] <- S
  }
  list(L = L_list[[r]], S = S_list[[r]], L_list = L_list, S_list = S_list, objective = objective)
}

#' Hard-thresholding operator
#'
#' @description
#' `hard_threshold()` implements the hard-thresholding operator on a given
#' matrix `D`, making `D` sparser: elements of `D` whose absolute value are less
#' than a given threshold `thresh` are set to 0, i.e. \eqn{D[|D| < thresh] = 0}.
#'
#' This is used in the non-convex PCP function [rrmc()] to provide a non-convex
#' replacement for the [prox_l1()] method used in the convex PCP function
#' [root_pcp()]. It is used to iteratively model the sparse `S` matrix with the
#' help of an adaptive threshold (`thresh` changes over the course of
#' optimization).
#'
#' @param D The input data matrix.
#' @param thresh The scalar-valued hard-threshold acting on `D` such that
#'   `D[i, j] = 0` when `abs(D[i, j]) < thresh`, and
#'   `D[i, j] = D[i, j]` otherwise.
#'
#' @returns The hard-thresholded matrix.
#'
#' @seealso [prox_l1()]
#' @examples
#' set.seed(42)
#' D <- matrix(rnorm(25), 5, 5)
#' S <- hard_threshold(D, thresh = 1)
#' D
#' S
#' @export
hard_threshold <- function(D, thresh) {
  D[abs(D) < thresh] <- 0
  D
}

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
