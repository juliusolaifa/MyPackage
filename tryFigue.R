vcov.MixMod <- function (object, parm = c("all", "fixed-effects", "var-cov","extra", 
                                          "zero_part"), sandwich = FALSE, ...) {
    parm <- match.arg(parm)
    V <- solve(object$Hessian)
    if (sandwich) {
        meat <- object$score_vect_contributions
        ind <- !names(meat) %in% "score.D" & !sapply(meat, is.null)
        meat[ind] <- lapply(meat[ind], rowsum, group = object$id[[1]], reorder = FALSE)
        meat <- do.call('cbind', meat)
        meat <- Reduce("+", lapply(split(meat, row(meat)), function (x) x %o% x))
        V <- V %*% meat %*% V
    }
    if (parm == "all") {
        return(V)
    }
    if (parm == "fixed-effects") {
        n_betas <- length(object$coefficients)
        return(V[seq_len(n_betas), seq_len(n_betas), drop = FALSE])
    }
    if (parm == "var-cov") {
        D <- object$D
        diag_D <- ncol(D) > 1 && all(abs(D[lower.tri(D)]) < sqrt(.Machine$double.eps))
        include <- if (diag_D) {
            unconstr_D <- log(diag(D))
            n_betas <- length(object$coefficients)
            seq(n_betas + 1, n_betas + length(unconstr_D))
        } else {
            unconstr_D <- chol_transf(D)
            n_betas <- length(object$coefficients)
            seq(n_betas + 1, n_betas + length(unconstr_D))
        }
        return(V[include, include, drop = FALSE])
    }
    if (parm == "extra") {
        if (is.null(object$phis)) {
            stop("the model behind 'object' contains no extra (phis) parameters.\n")
        } else {
            ind_phis <- grep("phi_", colnames(V), fixed = TRUE)
            return(V[ind_phis, ind_phis, drop = FALSE])
        }
    }
    if (parm == "zero_part") {
        if (is.null(object$gammas)) {
            stop("the fitted model does not have an extra zero part.")
        } else {
            gammas <- object$gammas
            ind_gammas <- grep("zi_", colnames(V), fixed = TRUE)
            return(V[ind_gammas, ind_gammas, drop = FALSE])
        }
    }
}
