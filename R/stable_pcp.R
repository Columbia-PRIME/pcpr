#' stablePCP
#'
#' \code{stablePCP} includes a non-negativity constraint on the \code{L} solution matrix.
#' It does not include LOD penalties. \cr \cr
#' Solve the following ADMM splitting problem: \cr \cr
#' min(L1,L2,L3,S1,S2) \cr
#' ||L1||_* + lambda * ||S1||_1 + mu/2 * ||L2+S2-D||_F^2 + I(L3>=0) \cr \cr
#' s.t. L1 = L2; L1 = L3; S1 = S2. \cr \cr
#'
#' @param D The original dataset.
#' @param lambda The \code{lambda} parameter penalizes the proximal L1 gradient on the \code{S} matrix.
#' @param mu The \code{mu} parameter penalizes the error term.
#'
#' @return Returns two solution matrices, the low rank \code{L} matrix and the sparse \code{S} matrix.
#'
#' @export
#'
stable_pcp <- function(D, lambda, mu) {
  pcp_lod(D, lambda, mu, LOD = 0)
}
