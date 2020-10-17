# is same function for convergence criteria
# is the difference among matrices > noise threshold?

# Compares L1, L2, L3 OR S1, S2
# Use ... for function to handle different number of inputs
# length(varargin) gives the number of function input arguments given in the call
is_same <- function(SAME_THRESH, ...) {
  flag <- TRUE
  varargin <- list(...)

  if (length(varargin) == 2) {
    if (max(abs(varargin[[1]] - varargin[[2]])) > SAME_THRESH) {
      flag <- FALSE
    }
  }
  else if (length(varargin) == 3) {
    if ((max(abs(varargin[[1]] - varargin[[2]])) > SAME_THRESH) |
        (max(abs(varargin[[1]] - varargin[[3]])) > SAME_THRESH) |
        (max(abs(varargin[[2]] - varargin[[3]])) > SAME_THRESH)) {
      flag <- FALSE
    }
  }
  flag
}
