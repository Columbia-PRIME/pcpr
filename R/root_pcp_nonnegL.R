#' Non-negative squareroot PCP function
#'
#' @export

#####
# Root PCP w/ non-negativity constraint
root_pcp_nonnegL <- function(D, lambda, mu) {

  n = nrow(D)
  p = ncol(D)
  rho = 0.1; # Augmented Lagrangian parameter

  # [L1,L2,L3,S1,S2,Z,Y1,Y2,Y3,Y4]
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

  MAX_ITER = 10000
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
    nuc = prox_nuclear( L2-Y1/rho, 1/rho )
    L1 = nuc[[1]]

    S1 = prox_l1( S2-Y2/rho, lambda/rho )
    Z = prox_fro( D-L2-S2-Y3/rho, mu/rho )
    L3 = pmax(L2-Y4/rho, 0)

    #% Update 2nd primal variable (L2,S2)
    term1 = L1+L3+D-Z + (Y1+Y4-Y3)/rho
    term2 = S1+D-Z + (Y2-Y3)/rho
    L2 = (2*term1 - term2) / 5
    S2 = (-term1 + 3*term2) / 5

    #% Update dual variable (Y1,Y2,Y3)
    Y1 = Y1 + rho*(L1-L2)
    Y2 = Y2 + rho*(S1-S2)
    Y3 = Y3 + rho*(L2+S2+Z-D)
    Y4 = Y4 + rho*(L3-L2)

    #%  Calculate primal & dual residuals; Update rho
    res_primal = sqrt( norm(L1-L2,'F')^2 +
                         norm(S1-S2,'F')^2 +
                         norm(Z+L2+S2-D,'F')^2 +
                         norm(L3-L2,'F')^2 )

    res_dual = rho * sqrt( 2 * norm(L2-L2_old,'F')^2 +
                             norm(S2-S2_old,'F')^2 +
                             norm(L2-L2_old+S2-S2_old,'F')^2 )

    if (res_primal > 10 * res_dual) {
      rho = rho * 2
    } else if (res_dual > 10 * res_primal) {
      rho = rho / 2}

    # %     % Calculate loss
    # %     loss(i) = nuclearL1 + lambda*sum(sum(abs(S1))) + mu*norm(L2+S2-D,'F') ...
    # %     + sum(sum(Y1.*(L1-L2))) + sum(sum(Y2.*(S1-S2))) ...
    # %     + sum(sum(Y3.*(L2+S2+Z-D))) + sum(sum(Y4.*(L3-L2))) ...
    # %     + rho/2 * ( sum(sum((L1-L2).^2)) + sum(sum((S1-S2).^2))
    # %     + sum(sum(L2+S2+Z-D)) + sum(sum(L3-L2)) );

    #% Check stopping criteria
    thresh_primal = EPS_ABS * sqrt(4*n*p) + EPS_REL *
      max(sqrt( norm(L1,'F')^2 + norm(S1,'F')^2 + norm(Z,'F')^2 ) + norm(L3,'F')^2,
          sqrt( 2 * norm(L2,'F')^2 + norm(S2,'F')^2 + norm(L2+S2,'F')^2 ),
          norm(D,'F') )

    thresh_dual = EPS_ABS * sqrt(4*n*p) + EPS_REL *
      sqrt( norm(Y1,'F')^2 + norm(Y2,'F')^2 + norm(Y3,'F')^2 + norm(Y4,'F')^2 )

    if (res_primal < thresh_primal && res_dual < thresh_dual) {
      flag_converge = 1
      print(paste0('Converged in ', i,' iterations.'))
      break}

  }

  L = (L1+L2+L3) / 3
  S = (S1+S2) / 2

  if (flag_converge == 0) {print('Did not converge.')}

  return(list(L,S))
}
