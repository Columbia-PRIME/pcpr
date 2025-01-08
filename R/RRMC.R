#' Given the observed mixture data matrix $D$, desired rank $r$,
#' and regularization hyperparameter $\eta$, RRMC seeks to find the
#' best estimates $\Lhat$ and $\Shat$ using an incremental rank-based
#' PCP strategy. That is, RRMC will first estimate a rank 1 model
#' $(\Lhat^{(1)}, \Shat^{(1)})$, before using that as the initialization
#' point to construct a rank 2 model $(\Lhat^{(2)}, \Shat^{(2)})$, and so on,
#' until the desired rank $r$ model $(\Lhat = \Lhat^{(r)}, \Shat = \Shat^{(r)})$
#' is recovered. All models from ranks 1 through $r$ are returned by RRMC in this way.
#' This incremental rank-based scheme is achieved via the following objective function:
#' \begin{equation} \label{eqn:rrmc}
#' \text{\rrmc}: \min_{\mb L,\mb S} \; \mathcal{I}_{\rank(\mb L) \leq r} + \eta \left\| \mb S \right\|_0 + \left\| \mb L + \mb S - \mb D \right\|_F^2.
#' \end{equation}
#' @export
RRMC <- function(D, r, eta, LOD = -Inf) {

  # 1. Initialize variables:
  # 1a. Matrix dimensions:
  n <- nrow(D)
  p <- ncol(D)

  # 1b. Process LOD input (convert LOD vector to matrix if needed, so it multiplies correctly):
  if (any(class(LOD) == "list")) LOD <- unlist(LOD)
  if (is.vector(LOD) && length(LOD) != 1) LOD <- matrix(LOD, n, p, T)

  # 1c. Omega = mask of observed values, for handling missingness later:
  Omega <- !is.na(D)
  D[!Omega] <- 0

  # 1d. mask_above/below_LOD = entries that are observed and above/below LOD:
  mask_above_LOD <- Omega & (D >= LOD)
  mask_below_LOD <- Omega & (D < LOD)

  # 1e. Hard-coded optimization parameters:
  epsilon <- 1e-05*norm(D, "F")
  temp <- svd(D)$d[1]
  t <- 10*log(40*r*n*temp/epsilon)
  zeta <- eta*temp

  # 1f. Instantiate variables to return to user:
  L <- matrix(0, nrow = n, ncol = p)
  L_list <- list()
  S_list <- list()
  objective <- vector("numeric", r*t)
  objective_index <- 1

  # 2. Iterative optimization procedure:
  for (k in 1:r) {
    for (i in 0:t) {
      S <- hard_thresholding(Omega*(D - L), zeta)
      D_hat <- L + S # approximation to the data
      # calculate the gradient of the LOD penalty:
      grad_LOD_penalty <- mask_above_LOD * (D_hat - D) + # observed and above LOD
                          mask_below_LOD * (D_hat < 0) * D_hat + # observed, below LOD, less than 0
                          ifelse(mask_below_LOD, (D_hat > LOD) * (D_hat - LOD), 0) # observed, below LOD, current appx bigger than LOD
                          #mask_below_LOD * (D_hat > LOD) * (D_hat - LOD) # observed, below LOD, current appx bigger than LOD

      # gradient step:
      M_i <- L - (n*p/sum(Omega))*grad_LOD_penalty
      L <- proj_rank_r(M_i, k)
      temp <- svd(M_i)$d
      zeta <- eta*(temp[k + 1] + 0.5^(i - 2)*temp[k])

      # calculate objective function value:
      objective[objective_index] <- norm(D - L - S, "F")^2 + eta*sum(S != 0)
      objective_index <- objective_index + 1
    }
    L_list[[k]] <- L
    S_list[[k]] <- S
  }

  list(L = L_list[[r]], S = S_list[[r]], L_list = L_list, S_list = S_list, objective = objective)

}

#' Hard-thresholding
#'
#' \code{hard_thresholding} implements hard-thresholding for the sparse matrix.
#'
#' @param S The input sparse matrix.
#' @param c The amount (or threshold) by which hard-thresholding penalizes \code{S}.
#'
#' @return The hard-thresholded sparse matrix.
#'
#' @seealso \code{\link{prox_l1}}
#' @examples
#' @keywords internal
hard_thresholding <- function(S, c) {

  X <- S
  X[abs(X) < c] <- 0
  X

}

#' Project rank
#'
#' \code{proj_rank_r} implements a best (i.e. closest) rank-r approximation of an input matrix.
#' This provides the non-convex replacement for nuclear norm on the \code{L} matrix, and is computed
#' by retaining the first r leading singular values/vectors of \code{L}. This is equivalent to
#' solving the following optimization problem: min || X - L ||_F s.t. rank(X) <= r.
#'
#' @param L The input low-rank matrix.
#' @param r The rank that \code{L} should be projected onto / truncated to.
#'
#' @return The best rank r approximation to the input matrix.
#'
#' @seealso \code{\link{prox_nuclear}}
#' @examples
#' @keywords internal
proj_rank_r <- function(L, r) {

  if (ncol(L) < r || nrow(L) < r) stop("r > matrix rank")
  if (ncol(L) == r || nrow(L) == r) return(L)

  USV <- svd(L)

  s <- USV$d
  s[(r+1):length(s)] <- 0
  S_new <- diag(s)

  X <- USV$u %*% S_new %*% t(USV$v)

  X

}
