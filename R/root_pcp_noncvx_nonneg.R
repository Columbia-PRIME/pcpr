#' Squareroot PCP function with non-convex replacement for nuclear norm and nonnegativity constraint
#'
#' \code{root_pcp_noncvx_nonneg} implements \code{rootPCP} with a non-negativity constraint on the \code{L} solution matrix and replaces the nuclear norm with a projection onto a smaller rank. \cr \cr
#' It solved the following ADMM splitting problem: \cr \cr
#' min_{L,S} \cr
#' 1_{rank(L) <= r}  + lambda * ||S||_1 + mu * ||L+S-D||_F \cr \cr
#' s.t. L >= 0. \cr \cr
#' This is first transformed to the problem: \cr \cr
#' min_{L1,L2,L3,S1,S2,Z} \cr
#'
#' 1_{rank(L) <= r} + lambda * ||S1||_1 + mu * ||Z||_F + I(L3>=0) \cr \cr
#' s.t. L1 = L2; S1 = S2; L2 + S2 + Z = D; L3 = L2 \cr \cr
#' The algorithm conducts ADMM splitting as (L1,S1,Z,L3),(L2,S2).
#'
#' @param D The original dataset.
#' @param lambda The \code{lambda} parameter penalizes the proximal L1 gradient on the \code{S} matrix.
#' @param mu The \code{mu} parameter penalizes the error term.
#' @param r The \code{r} parameter specifies the desired rank.
#' @param verbose A logical indicating if you would like information on the number of iterations required to reach convergence printed. Optional, and by default \code{verbose = FALSE}.
#'
#' @return Returns two solution matrices, the low rank \code{L} matrix and the sparse \code{S} matrix.
#'
#' @export
root_pcp_noncvx_nonneg <- function(D, lambda, mu, r, verbose=FALSE) {

  n = nrow(D)
  p = ncol(D)
  rho = 0.1; # Augmented Lagrangian parameter

  L1 <- matrix(0, n, p)
  L2 <- matrix(0, n, p)
  L3 <- matrix(0, n, p)

  S1 <- matrix(0, n, p)
  S2 <- matrix(0, n, p)

  Z  <- matrix(0, n, p)
  Y1 <- matrix(0, n, p)
  Y2 <- matrix(0, n, p)
  Y3 <- matrix(0, n, p)
  Y4 <- matrix(0, n, p)

  MAX_ITER = 20000
  EPS_ABS = 1e-6
  EPS_REL = 1e-6

  flag_converge = 0
  #% loss = zeros(MAX_ITER, 1);

  #% ADMM-splitting iterations
  for (i in 1:MAX_ITER) {

    #% Store previous values of L2,S2
    L2_old = L2
    S2_old = S2

    #% Update 1st primal variable (L1,S1,Z,L3)
    ## This makes it non-convex
    L1 = proj_rank_r( L2-Y1/rho, r)

    S1 = prox_l1( S2-Y2/rho, lambda/rho )
    Z = prox_fro( D-L2-S2-Y3/rho, mu/rho )
    L3 = pmax(L2-Y4/rho, 0)

    #% Update 2nd primal variable (L2,S2)
    L2 = 1/3*( D-Z+ 2*L1 -S1 + (2*Y1 -Y2-Y3)/rho )
    S2 = 1/3*( D-Z+ 2*S1 -L1 + (2*Y2 -Y1-Y3)/rho )

    #% Update dual variable (Y1,Y2,Y3)
    Y1 = Y1 + rho*(L1-L2)
    Y2 = Y2 + rho*(S1-S2)
    Y3 = Y3 + rho*(L2+S2+Z-D)

    #%  Calculate primal & dual residuals; Update rho
    res_primal = sqrt( norm(L1-L2,'F')^2 +
                         norm(S1-S2,'F')^2 +
                         norm(Z+L2+S2-D,'F')^2 +
                         norm(L3-L2,'F')^2 )
    res_dual = rho * sqrt(2*norm(L2-L2_old,'F')^2 +
                             norm(S2-S2_old,'F')^2 +
                             norm(L2-L2_old+S2-S2_old,'F')^2 )

    if (res_primal > 10 * res_dual) {
      rho = rho * 2
    } else if (res_dual > 10 * res_primal) {
      rho = rho / 2}

    # %     % Calculate loss
    # %     loss(i) = nuclearL1 + lambda*sum(sum(abs(S1))) + mu*norm(L2+S2-D,'fro') ...
    # %         + sum(sum(Y1.*(L1-L2))) + sum(sum(Y2.*(S1-S2))) + sum(sum(Y3.*(L2+S2+Z-D))) + ...
    # %          sum(sum(Y4.*(L3-L2))) + rho/2 * ( sum(sum((L1-L2).^2)) + sum(sum((S1-S2).^2)) + ...
    # %          sum(sum(L2+S2+Z-D)) + sum(sum(L3-L2)) );

    #% Check stopping criteria
    thresh_primal = EPS_ABS * sqrt(4*n*p) + EPS_REL * max(
      sqrt( norm(L1,'F')^2 + norm(S1,'F')^2 + norm(Z,'F')^2) + norm(L3,'F')^2,
      sqrt( 2*norm(L2,'F')^2 + norm(S2,'F')^2 + norm(L2+S2,'F')^2 ),
      norm(D,'F'))

    thresh_dual = EPS_ABS * sqrt(4*n*p) + EPS_REL *
      sqrt( norm(Y1,'F')^2 + norm(Y2,'F')^2 + norm(Y3,'F')^2  + norm(Y4,'F')^2 )

    if (res_primal < thresh_primal && res_dual < thresh_dual) {
      flag_converge = 1
      if (verbose) print(paste0('Converged in ', i,' iterations.'))
      break}

  }

  L = pmax((L1+L2+L3) / 3, 0)
  S = (S1+S2) / 2

  if (flag_converge == 0 & verbose) print('Did not converge.')

  return(list(L=L,S=S))
}

