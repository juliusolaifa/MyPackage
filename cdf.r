bivarnormcdf <- function(sigma) {
    sig11 <- sqrt(sigma[1,1])
    sig22 <- sqrt(sigma[2,2])
    sig12 <- sigma[1,2]
    rho <- sig12 / (sig11 * sig22)
    a <- sig12/sig11^2
    if(rho > 0) {
        lower = -Inf
        upper = 0
    } else {
        lower = 0
        upper = Inf
    }
    bivarnormpdf <- function(x, y, sig11, sig22, rho) {
      exponent <- -(1 / (2 * (1 - rho^2))) * (x^2 / sig11^2 + y^2 / sig22^2 - 2 * rho * x * y / (sig11 * sig22))
      density <- (1 / (2 * pi * sig11 * sig22 * sqrt(1 - rho^2))) * exp(exponent)
      return(density)
    }
    result <- integrate(function(y) {
        sapply(x, function(y_val) {
            integrate(function(x) bivarnormpdf(x, y_val, sig11, sig22, rho), lower = y_val/a, upper = 0)$value
        })
    }, lower = lower, upper = upper)
    # Check for errors
    if (result$message != "OK") {
        warning("Integration did not converge")
    }
    return(result$value)
}
mat <- matrix(c(2,1,1,2),2)
bivarnormcdf(mat)
