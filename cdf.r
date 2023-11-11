bivarnormcdf <- function(sigma, lower, upper) {
    sig11 <- sqrt(sigma[1,1])
    sig22 <- sqrt(sigma[2,2])
    sig12 <- sigma[1,2]
    rho <- sig12 / (sig11 * sig22)

    bivarnormpdf <- function(x, y, sig11, sig22, rho) {
      exponent <- -(1 / (2 * (1 - rho^2))) * (x^2 / sig11^2 + y^2 / sig22^2 - 2 * rho * x * y / (sig11 * sig22))
      density <- (1 / (2 * pi * sig11 * sig22 * sqrt(1 - rho^2))) * exp(exponent)
      return(density)
    }
    
    result <- integrate(function(x) {
        sapply(x, function(x_val) {
            integrate(function(y) bivarnormpdf(x_val, y, sig11, sig22, rho), lower = lower[1], upper = upper[1])$value
        })
    }, lower = lower[2], upper = upper[2])

    # Check for errors
    if (result$message != "OK") {
        warning("Integration did not converge")
    }
    
    return(result$value)
}


bivarnormcdf(matrix(c(2,1,1,2),2))