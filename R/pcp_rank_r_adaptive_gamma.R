#' @export
pcp_rank_r_adaptive_gamma = function(D,gamma,r, L_init = NULL, verbose=FALSE) {

    # % solve the problem
    # %
    # %   min I_{rank(L) <=r} + gamma ||S||_1 + .5*||PO(Y-L-S)||_F^2
    # %
    # % using a proximal gradient method

    n = nrow(D)
    p = ncol(D)

    Om = !is.na(D)
    D[!Om] = 0

    L = matrix(0, nrow=n, ncol=p)
    S = matrix(0, nrow=n, ncol=p)

    if (is.null(L_init)) {
        L = (n*p / sum(Om)) * proj_rank_r(D,r) # best rank r approximation as initialization point.
    } else {
        L = (n*p / sum(Om)) * proj_rank_r(L_init,r) # best rank r approximation as initialization point.
    }

    obj = Inf

    t = .01
    iter = 0
    done = FALSE
    max_iter = 4000
    allObj = c()
    gamma_vec = apply(Om, 1, function(v) {
        denom = sum(v) - r
        if (denom <= 0) denom = 1
        sqrt((n - r) / denom)
    })

    while (!done) {

        R = Om * ( D - L - S )
        S_new = prox_l1( S + t * R, t*gamma*gamma_vec )
        L_new = proj_rank_r( L + t * R, r )

        delta = sqrt(norm(S-S_new,'F')^2 + norm(L-L_new,'F')^2)

        # norm 1 is the maximum absolute column sum of the matrix.
        obj_new = gamma * norm(S,"1") + .5 * norm( Om * ( D - L - S ), 'F' )^2

        if (obj_new > obj) {
            t = t * .95
        } else {
            t = t * 1.01
            S = S_new
            L = L_new
            delta_obj = obj - obj_new;
            obj = obj_new}

        iter = iter + 1
        allObj = c(allObj, obj)

        if (verbose & iter %% 100 == 0) {print(paste0('Iter: ', iter, ', Obj: ', round(obj,4), ', Step: ', round(delta,4)))}

        if (iter >= max_iter | t < 1e-10) {done = TRUE}
        if (iter >= max_iter & verbose) {print('Did not converge.')}
        if (verbose & t < 1e-10) {print(paste("Converged in:", iter, "iterations."))}
    }
        return(list(L = L, S = S, final_iter = iter))
}





