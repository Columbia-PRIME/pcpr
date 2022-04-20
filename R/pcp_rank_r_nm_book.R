#' @export
pcp_rank_r_nm_book = function(D, gamma, r, verbose = FALSE) {
  #% solve the problem
  #%
  #%   min I_{rank(L) <=r} + gamma ||S||_1 + .5*||PO(Y-L-S)||_F^2
  #%
  #% using a nesterov's accelerated proximal gradient method
  n = nrow(D)
  p = ncol(D)
  Om = !is.na(D)
  D[!Om] = 0

  L = matrix(0, nrow=n, ncol=p)
  S = matrix(0, nrow=n, ncol=p)

  L = (n*p / sum(Om)) * proj_rank_r(D,r) # best rank r approximation as initialization point.

  # set up for nesterov's method:
  x_L = L
  x_S = S
  v_L = L
  v_S = S
  x_L_prev = L
  x_S_prev = S

  obj = Inf
  t = 1
  iter = 0
  done = FALSE
  max_iter = 10000
  allObj = c()

  while (!done) {
    
    t_new = (1 + sqrt(1 + 4 * t^2)) / 2
    beta = t / t_new
    
    x_L_prev = x_L
    x_S_prev = x_S

    # apply a prox gradient step at y
    R = Om * ( D - v_L - v_S )
    x_S = prox_l1(v_S + t * R, t*gamma)
    x_L = proj_rank_r(v_L + t * R, r)

    v_L = x_L + beta * (x_L - x_L_prev)
    v_S = x_S + beta * (x_S - x_S_prev)

    obj_x = gamma * sum(abs(x_S)) + 0.5 * norm(Om * (D - x_L - x_S), 'F')^2
    #obj_v = gamma * sum(abs(v_S)) + 0.5 * norm(Om * (D - v_L - v_S), 'F')^2
    
    obj = obj_x

    iter = iter + 1
    t = t_new

    delta = sqrt(norm(x_L - x_L_prev, 'F')^2 + norm(x_S + x_S_prev, 'F')^2)

    if (iter >= max_iter | t < 1e-10) done = TRUE
    if (verbose) {
      if (iter %% 100 == 0) print(paste0('Iter: ', iter, ', Obj: ', round(obj,4), ', Step: ', round(delta,4)))
      if (iter >= max_iter) print('Did not converge.')
      if (t < 1e-10 & iter < max_iter) print(paste("Converged in:", iter, "iterations."))
    }

    allObj = c(allObj, obj)
  }
  return(list(L = x_L, S = x_S, final_iter = iter, obj_vals = allObj))
}