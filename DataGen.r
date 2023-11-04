herit.glmmdata <- function(num.iterations, x, beta, sigma.u, cluster.sizes, family, ...) {
  
  # E(y|u) = log(xb + zu), u ~ N(0,E)
  cond_mean <- function() {
    if(!is.matrix(sigma.u)) sigma.u <- as.matrix(sigma.u)
    x <- model.matrix(~x)
    u <- mvrnorm(n=length(cluster.sizes), mu=rep(0,ncol(sigma.u)), Sigma = sigma.u)
    group <- factor(rep(1:length(cluster.sizes), cluster.sizes))
    z <- model.matrix(~ 0 + group)
    z <- do.call(cbind, lapply(1:ncol(u), function(i) z * x[, i]))
    return(log(x%*%beta + z%*%c(u)))
  }

  total_observations <- sum(cluster.sizes)
  y_matrix <- matrix(nrow=num.iterations, ncol=total_observations)

  generate <- switch(family,
           nbinom2=MASS::rnegbin,
           tweedie=tweedie::rtweedie,
           compois=COMPoissonReg::rcmp,
           genpois=HMMpa::rgenpois,
           stop("data generation not implemented for family: ", family)
  )
                               
  for (i in 1:num.iterations) {                            
      #split conditional means into list according to group
      cond_means <- split(cond_mean(), rep(1:length(cluster.sizes), cluster.sizes))
      y <- unlist(mapply(generate, n=cluster.sizes, mu=cond_means, ...))
      y_matrix[i, ] <- y
  }

  result <- structure(list("x" = x, "y" = y_matrix, "family" = family), class = "heritData")
  return(result)
}
                               

# Example
x <- c(23,24,25,26,29,30)
beta <- c(5,9)
sigma.u <- 3#matrix(c(2,1,1,2),2)
cluster.sizes <- c(3,2,1)
n = 2
family = "nbinom2"
theta = 2
num.iterations =1
herit.glmmdata(num.iterations, x, beta, sigma.u, cluster.sizes, family, theta)
