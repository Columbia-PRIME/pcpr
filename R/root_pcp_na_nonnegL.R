#' Nonnegative squareroot PCP function with missing values (NA)
#'
#' \code{root_pcp_na_nonnegL} implements \code{rootPCP} with a non-negativity constraint on the \code{L} solution matrix. \cr \cr
#' It solved the following ADMM splitting problem: \cr \cr
#' min(L,S) \cr
#' ||L||_* + lambda * ||S||_1 + mu * || P_(obs)[L+S-D] ||_F + I[L>=0] \cr \cr
#' This is first transformed to the problem: \cr \cr
#' min(L1,L2,L3,S1,S2,Z) \cr
#' ||L1||_* + lambda * ||S1||_1 + mu * ||Z||_F + I[L3>=0] \cr \cr
#' s.t. L1 = L2; S1 = S2; Z = P_obs[D - L2 - S2]; L1 = L3. \cr \cr
#' The algorithm conducts ADMM splitting as (L1,S1,Z),L3,(L2,S2). \cr \cr
#' This version allows for missing values. \cr \cr
#' Use NA for missing entries in D. \cr \cr
#' Assume that the true L>=0
#'
#' @param D The original dataset.
#' @param lambda The \code{lambda} parameter penalizes the proximal L1 gradient on the \code{S} matrix.
#' @param mu The \code{mu} parameter penalizes the error term.
#' @param verbose A logical indicating if you would like information on the number of iterations required to reach convergence printed. Optional, and by default \code{verbose = FALSE}.
#'
#' @return Returns two solution matrices, the low rank \code{L} matrix and the sparse \code{S} matrix.
#'
#' @export
root_pcp_na_nonnegL <- function(D, lambda, mu, verbose = FALSE) {

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

# mask: support of observation of D
mask = !is.na(D)
D[!mask] = 0

MAX_ITER = 10000
EPS_ABS = 1e-6
EPS_REL = 1e-6

flag_converge = 0
#% loss = zeros(MAX_ITER, 1);

#% ADMM-splitting iterations
for (i in 1:MAX_ITER) {

#% Store previous values of L2,S2
L2_old = L2
L3_old = L3
S2_old = S2

#% Update 1st primal variable (L1,S1,Z)
nuc = prox_nuclear( (L2+L3-Y1/rho-Y4/rho)/2, 1/rho/2  )
L1 = nuc[[1]]

S1 = prox_l1( S2-Y2/rho, lambda/rho )
Z = prox_fro( mask*(D-L2-S2)-Y3/rho, mu/rho )
L3 = pmax( L1 + Y4/rho, 0)

#% Update 2nd primal variable (L2,S2)
L2_obs = mask * (1/3 *( D - Z + 2*L1 - S1 + (2*Y1 - Y2 - Y3) / rho ))
L2_unobs = (1-mask) * (L1+Y1/rho)
L2 = L2_obs + L2_unobs

S2_obs = mask * (1/3 * ( D - Z + 2*S1 - L1 + (2*Y2 - Y1 - Y3) / rho ))
S2_unobs = (1-mask) * (S1+Y2/rho)
S2 = S2_obs + S2_unobs

#% Update dual variable (Y1,Y2,Y3)
Y1 = Y1 + rho*(L1-L2)
Y2 = Y2 + rho*(S1-S2)

Y3 = Y3 + rho * (Z - mask*(D-L2-S2))
Y4 = Y4 + rho * (L1 - L3)

#%  Calculate primal & dual residuals; Update rho
res_primal = sqrt( norm(L1-L2,'F')^2 +
                     norm(S1-S2,'F')^2 +
                       norm(Z-mask*(D-L2-S2),'F')^2 + norm(L1-L3,'F')^2)

res_dual = rho * sqrt( norm(L2+L3-L2_old-L3_old,'F')^2 +
                         norm(S2-S2_old,'F')^2 +
                           norm(mask * (L2-L2_old+S2-S2_old),'F')^2 )

if (res_primal > 10 * res_dual) {
  rho = rho * 2
  } else if (res_dual > 10 * res_primal) {
    rho = rho / 2}

#% Check stopping criteria
thresh_primal = EPS_ABS * sqrt(4*n*p) + EPS_REL *
                max(sqrt( 2*norm(L1,'F')^2 + norm(S1,'F')^2 + norm(Z,'F')^2 ),
                    sqrt( norm(L2,'F')^2 + norm(S2,'F')^2 + norm(mask*(L2+S2),'F')^2 + norm(L3,'F')^2),
                    norm(D,'F'))

thresh_dual = EPS_ABS * sqrt(3*n*p) + EPS_REL *
              sqrt( norm(Y1+Y4,'F')^2 + norm(Y2,'F')^2 + norm(Y3,'F')^2)

if (res_primal < thresh_primal && res_dual < thresh_dual) {
  flag_converge = 1
  if (verbose) print(paste0('Converged in ', i,' iterations.'))
  break}

}

L = (L1+L2+L3) / 3
S = (S1+S2) / 2

if (flag_converge == 0 & verbose) print('Did not converge.')

return(list(L=L,S=S))
}

