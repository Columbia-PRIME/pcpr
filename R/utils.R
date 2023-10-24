#' Retrieves the default PCP parameter settings for a given matrix.
#'
#' \code{get_pcp_defaults} retrieves the default settings of \code{lambda}, \code{mu}, and \code{eta} for a given matrix \code{D}.
#'
#' @param D The matrix for which the default \code{lambda}, \code{mu}, and \code{eta} are needed.
#'
#' @return A list containing the default lambda and mu values used in \code{\link{root_pcp}}, and the default eta value used in \code{\link{RRMC}}. Labelled as "lambda", "mu", and "eta" respectively.
#' @examples
#'
#' # Examine the queens PM2.5 data
#' queens
#'
#' # Get rid of the Date column
#' D <- queens[, 2:ncol(queens)]
#'
#' # Get the default PCP parameters
#' get_pcp_defaults(D)
#' @export
get_pcp_defaults <- function(D) {
  n <- nrow(D)
  p <- ncol(D)
  l <- 1/sqrt(max(n, p))
  m <- sqrt(min(n, p)/2)
  e <- l/m
  list(lambda = l, mu = m, eta = e)
}
