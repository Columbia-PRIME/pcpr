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