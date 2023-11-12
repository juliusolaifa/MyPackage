bivarnormcdf <- function(sigma, lower="-Inf", upper=Inf, inner_lower="-Inf", inner_upper=Inf) {

    bivarnormpdf <- function(x1, x2, sig11, sig22, rho) {
        exponent <- -(1 / (2 * (1 - rho^2))) * (x1^2 / sig11^2 + x2^2 / sig22^2 - 2 * rho * x1 * x2 / (sig11 * sig22))
        density <- (1 / (2 * pi * sig11 * sig22 * sqrt(1 - rho^2))) * exp(exponent)
        return(density)
    }

    if(!is.character(lower) && is.character(inner_lower)) { #x2 depends on x1
        sigma <- sigma
    }else if(is.character(lower) && !is.character(inner_lower)) { #x depends on x2
        #swap the diagonal of the matrix and swap the integral
        temp <- sigma[1,1]; sigma[1,1] <- sigma[2,2]; sigma[2,2] <- temp
        temp2 <- lower; lower <- inner_lower; inner_lower <- temp2
        temp3 <- upper; upper <- inner_upper; inner_upper <- temp3
    } else if(is.character(lower) && is.character(inner_lower)) {
        stop("Only one of `lower` and `inner_lower can be passed as string")
    }
    sig11 <- sqrt(sigma[1,1]); sig22 <- sqrt(sigma[2,2]); sig12 <- sigma[1,2]
    rho <- sig12 / (sig11 * sig22)
    a <- sig12 / sig11^2
  
    result <- integrate(function(x1) {
        env <- environment()
        sapply(x1, function(x1_val) {
            # Evaluate the string expression for inner_lower
            eval_inner_lower <- eval(parse(text=gsub("x1", as.character(x1_val), inner_lower)), envir=env)
            integrate(function(x2) bivarnormpdf(x1_val, x2, sig11, sig22, rho), 
                      lower = eval_inner_lower, 
                      upper = inner_upper)$value
        })
    }, lower = lower, upper = upper)

    # Check for errors
    if (result$message != "OK") {
        warning("Integration did not converge")
    }
    return(result$value)
}
mat <- matrix(c(4,2,2,3),2)
bivarnormcdf(mat, lower=0, upper=Inf, inner_lower=0, inner_upper=Inf)
#bivarnormcdf(mat, lower="a*x1", upper=Inf, inner_lower=-Inf, inner_upper=0)
