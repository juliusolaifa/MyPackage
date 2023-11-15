vcov2 <- function(modObj) {
    beta <- unname(modObj$coefficients)
    sigma_chol <- GLMMadaptive:::chol_transf(modObj$D)
    phis <- unname(modObj$phis)
    theta_chol <- c(beta, sigma_chol, phis)

    nbeta <- length(beta)
    nsigma <- length(sigma_chol)
    sigma_f_ind <- nbeta + 1
    sigma_l_ind <- nbeta + nsigma 
    
    V_chol <- vcov(modObj)
    V <- delta(theta_chol, V_chol, function(theta) {
        sigma <- theta[sigma_f_ind:sigma_l_ind]
        sigma <- GLMMadaptive:::chol_transf(sigma)
        sigma <- sigma[upper.tri(sigma, diag = TRUE)]
        theta[sigma_f_ind:sigma_l_ind] <- sigma
        theta
    }) 
    dimnames(V) <- dimnames(V_chol)
    theta <- chol_tocov(theta_chol)
    
    return(list("V" = V, "theta" = theta, "sigma_f_ind" = sigma_f_ind, "nsigma" = nsigma))
}

vcov.heritMod <- function(modObj) {
    V <- vcov2(modObj)
    theta <- V$theta
    nsigma <- V$nsigma
    sigma_f_ind <- V$sigma_f_ind
    V <- V$V
    
    
    if (nsigma == 1) {
        V1 <- V2 <- V
        V2[sigma_f_ind, ] <- V2[, sigma_f_ind] <- 0
        pi_1 <- 1; pi_2 <- 0.5
        V.combine <- pi_1V1 + pi_2*V2
    } else if(nsigma== 3) {
        pi_1 <- 0.25; pi_2 <- 0.25; pi_3 <- 0.25; pi_4 <- 0.25
        V1 <- V2 <- V3 <- V4 <- V
        V2[sigma_f_ind, ] <- V2[, sigma_f_ind] <- 0
        V3[c(sigma_f_ind+2), ] <- V3[, c(sigma_f_ind+2)] <- 0
        V4[c(sigma_f_ind,sigma_f_ind+2), ] <- V4[,c(sigma_f_ind,sigma_f_ind+2)] <- 0
        V.combine <- pi_1*V1 + pi_2*V2 + pi_3*V3 + pi_4*V4
    }
    delta(theta,V.combine,****) 
}

delta <- function(estimates, cov_matrix, transformation) {
    if (!is.vector(estimates) || !is.matrix(cov_matrix)) {
        stop("Invalid inputs: 'estimates' must be a vector and 'cov_matrix' must be a matrix.")
    }
    grad <- numDeriv::jacobian(func = transformation, x = estimates)
    transformed_var_cov <- grad %*% cov_matrix %*% t(grad)
    if (length(transformation(estimates)) == 1) {
        return(transformed_var_cov[1, 1])
    } else {
        return(transformed_var_cov)
    }
}

confint.heritMod <- function() {
}
