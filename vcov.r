vcov2 <- function(modObj) {
    beta <- unname(modObj$coefficients)
    sigma_chol <- GLMMadaptive:::chol_transf(modObj$D)
    phis <- unname(modObj$phis)
    theta_chol <- c(beta, sigma_chol, phis)

    sigma_f_ind <- length(beta) + 1
    sigma_l_ind <- length(beta) + length(sigma_chol)
    
    chol_tocov <- function(theta) {
        sigma <- theta[sigma_f_ind:sigma_l_ind]
        sigma <- GLMMadaptive:::chol_transf(sigma)
        sigma <- sigma[upper.tri(sigma, diag = TRUE)]
        theta[sigma_f_ind:sigma_l_ind] <- sigma
        theta
    }
    
    J <- numDeriv::jacobian(chol_tocov, theta_chol)
    V_chol <- vcov(modObj)
    V <- J %*% V_chol %*% t(J)
    dimnames(V) <- dimnames(V_chol)
    V
}

vcov.heritMod <- function(modObj) {
    V <- vcov2(modObj)
}

confint.heritMod <- function() {
}
