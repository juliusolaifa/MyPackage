

cond_mean <- function(x, beta, sigma.u, cluster.sizes) {
    if(!is.matrix(sigma.u)) sigma.u <- as.matrix(sigma.u)
    x <- model.matrix(~x)
    u <- mvrnorm(n=length(cluster.sizes), mu=rep(0,ncol(sigma.u)), Sigma = sigma.u)
    group <- factor(rep(1:length(cluster.sizes), cluster.sizes))
    z <- model.matrix(~ 0 + group)
    z <- do.call(cbind, lapply(1:ncol(u), function(i) z * x[, i]))
    return(log(x%*%beta + z%*%c(u)))
}

# Example
library(MASS)
x <- c(23,24,25,26,29,30)
beta <- c(5,9)
sigma.u <- 3#matrix(c(2,1,1,2),2)
cluster.sizes <- c(3,2,1)
cond_mean(x, beta, sigma.u, cluster.sizes)
