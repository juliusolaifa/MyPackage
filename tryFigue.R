vcov.MixMod <- function (object) {
    V <- solve(object$Hessian)
    return(V)
}

# central difference vector (or matrix) for a function
cd_vec <- function (x, f, ..., eps = 0.001) {
    n <- length(x)
    res <- matrix(0, n, n)
    ex <- pmax(abs(x), 1)
    for (i in seq_len(n)) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[, i] <- diff.f / diff.x
    }
    0.5 * (res + t(res))
}


#tht: These are the parameters of interest
list_thetas <- list(betas = betas, D = if (diag_D) log(diag(D)) else chol_transf(D))
    if (!is.null(phis)) {
        list_thetas <- c(list_thetas, list(phis = phis))
    }
               
tht <- unlist(as.relistable(list_thetas))

Hessian <- cd_vec(tht, score_mixed, id = id, y = y, N = N, X = X, Z = Z, 
                      offset = offset, X_zi = X_zi, Z_zi = Z_zi, offset_zi = offset_zi, 
                      GH = GH, canonical = canonical, user_defined = user_defined, 
                      Xty = Xty, Xty_weights = Xty_weights, log_dens = log_dens, mu_fun = mu_fun, 
                      var_fun = var_fun, mu.eta_fun = mu.eta_fun, 
                      score_eta_fun = score_eta_fun, score_eta_zi_fun = score_eta_zi_fun, 
                      score_phis_fun = score_phis_fun, list_thetas = list_thetas, 
                      diag_D = diag_D, penalized = penalized, pen_mu = pen_mu, 
                      pen_invSigma = pen_invSigma, pen_df = pen_df, weights = weights)
