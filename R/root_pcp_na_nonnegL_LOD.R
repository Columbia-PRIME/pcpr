#' Squareroot PCP function with missing values (NA)
#'
#' \code{root_pcp} implements \code{rootPCP} with NO non-negativity constraint on the \code{L} solution matrix. \cr \cr
#' It solved the following ADMM splitting problem: \cr \cr

#' min_{L,S} \cr
#' ||L||_* + lambda * ||S||_1 + mu * (\sum_{ij observed} f((L+S)_ij, D_ij))^0.5 + I{L>=0} \cr \cr
#' This is first transformed to the problem: \cr \cr
#'
#' min_{L1,L2,S,Z} \cr
#' ||L1||_* + lambda * ||S||_1 + mu * (\sum_{ij observed} f(Z_ij,D_ij))^0.5 + I{L2>=0} \cr \cr
#' s.t. L1 = L2; L1 + S = Z. \cr \cr
#' The algorithm conducts ADMM splitting as L1, S, Z, L2. \cr \cr
#' This version allows for missing values. \cr \cr
#' Use NA for missing entries in D. \cr \cr
#' Assume that the true L >= 0 & observations in D >= 0 \cr \cr
#' Use -1 for <LOD \cr \cr
#'
#' LOD penalty: \cr \cr
#' f(x,y) = (x-y)^2 if y>0 \cr
#'        = (x-Delta)^2 if y=-1, x>Delta \cr
#'        = x^2 if y=-1, x<0 \cr
#'        = 0 otherwise
#' @param D The original dataset.
#' @param lambda The \code{lambda} parameter penalizes the proximal L1 gradient on the \code{S} matrix.
#' @param mu The \code{mu} parameter penalizes the error term.
#'
#' @return Returns two solution matrices, the low rank \code{L} matrix and the sparse \code{S} matrix.
#'
#' @export
#'
root_pcp_na_nonnegL_LOD <- function(D, lambda, mu, Delta) {

n = nrow(D)
p = ncol(D)
rho = 0.1; # Augmented Lagrangian parameter

L1 <- matrix(0, n, p)
L2 <- matrix(0, n, p)

S <- matrix(0, n, p)

Z  <- matrix(0, n, p)
Y1 <- matrix(0, n, p)
Y2 <- matrix(0, n, p)

# mask: support of observation of D

mask_above_lod = D >= 0
mask_above_lod[is.na(mask_above_lod)] = 0

mask_below_lod = D < 0
mask_below_lod[is.na(mask_below_lod)] = 0

mask_obs = !is.na(D)
D[!mask_obs] = -2

MAX_ITER = 10000
EPS_ABS = 1e-6
EPS_REL = 1e-6

flag_converge = 0

#% ADMM-splitting iterations
for (i in 1:MAX_ITER) {

#% Store previous values of L2,S<Z
L2_old = L2
S_old = S
Z_old = Z

#% Update 1st primal variable (L1,S1,Z)
nuc = prox_nuclear( (L2+Z-S-Y1/rho-Y2/rho)/2, 1/rho/2 )
L1 = nuc[[1]]

S = prox_l1( Z-L1-Y2/rho, lambda/rho )

temp = L1+S+Y2/rho

Z_unobs = (1-mask_obs) * temp
Z_obs_below_LOD1 = (mask_below_lod & (temp>=0) & (temp<=Delta)) * temp

temp2 = mask_obs * (1-(mask_below_lod & (temp>=0) & (temp<=Delta))) * temp -
         (mask_above_lod * D) - (Delta * mask_below_lod * (temp>=Delta))

Z = prox_fro( temp2, mu/rho ) + (mask_above_lod * D) +
    (Delta * mask_below_lod * (temp>=Delta)) +
    Z_unobs + Z_obs_below_LOD1

L2 = pmax(L1+Y1/rho,0)

#% Update dual variable (Y1,Y2)
Y1 = Y1 + rho*(L1-L2)
Y2 = Y2 + rho*(L1+S-Z)

#%  Calculate primal & dual residuals; Update rho
res_primal = sqrt(norm(L1-L2,'F')^2 + norm(L1+S-Z,'F')^2)
res_dual   = rho * sqrt( norm(L2-L2_old+Z-Z_old-S+S_old,'F')^2 + norm(Z-Z_old,'F')^2)

if (res_primal > 10 * res_dual) {
  rho = rho * 2
  } else if (res_dual > 10 * res_primal) {
    rho = rho / 2}

#% Check stopping criteria
thresh_primal = EPS_ABS * sqrt(3*n*p) + EPS_REL *
                max( sqrt( norm(L1,'F')^2 + norm(L1+S,'F')^2),
                     sqrt( norm(L2,'F')^2 + norm(Z,'F')^2) )

thresh_dual = EPS_ABS * sqrt(3*n*p) + EPS_REL *
               sqrt( norm(Y1+Y2,'F')^2 + norm(Y2,'F')^2 )

if (res_primal < thresh_primal && res_dual < thresh_dual) {
  flag_converge = 1
  print(paste0('Converged in ', i,' iterations.'))
  break}

}

L_final = (L1+L2) / 2
S_final = S

if (flag_converge == 0) {print('Did not converge.')}

return(list(L = L_final, S = S_final))
}

