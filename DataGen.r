herit.glmmdata <- function(iter, x, beta, sigma.u, cluster.sizes, family, ...) {
  sigma.u <- as.matrix(sigma.u)
  x <- as.matrix(x)
  if(nrow(x) != sum(cluster.sizes)) stop("number of rows of x must be equal to sum of cluster sizes")
  x <- model.matrix(~x)
  u <- mvrnorm(n=length(cluster.sizes), mu=rep(0,ncol(sigma.u)), Sigma = sigma.u)
  group <- factor(rep(1:length(cluster.sizes), cluster.sizes))
  z0 <- model.matrix(~ 0 + group)
  z <- do.call(cbind, lapply(1:ncol(u), function(i) z0 * x[, i]))

  total_observations <- sum(cluster.sizes)
  y_matrix <- matrix(nrow=iter, ncol=total_observations)

  cluster_assignments <- rep(paste0("strain", 1:length(cluster.sizes)), times=cluster.sizes)
  cluster_assignments <- factor(cluster_assignments)
                              
  #split conditional means into list according to group
  cond_means <- split(log(x%*%beta + z%*%c(u)), rep(1:length(cluster.sizes), cluster.sizes))
                               
  args <- list(...) 
  generate <- switch(family,
           nbinom2={
             #MASS::rnegbin
             function() unlist(mapply(MASS::rnegbin, n=cluster.sizes, mu=cond_means,args$theta))
           },
           tweedie={
             #tweedie::rtweedie
             function() unlist(mapply(tweedie::rtweedie, n=cluster.sizes, mu=cond_means, args$phi, args$power))
           },
           compois={
             #COMPoissonReg::rcmp
             function() unlist(mapply(COMPoissonReg::rcmp, n=cluster.sizes, lambda=cond_means, args$nu))
           },
           # genpois={
           #   #HMMpa::rgenpois
           #   lambda1 <- unlist(sapply(cond_means, function(cm) cm *(1 - args$lambda2)))
           #   function() unlist(mapply(HMMpa::rgenpois, n=cluster.sizes, lambda1, args$lambda2))
           # },
           stop("data generation not implemented for family: ", family)
  )

                            
  for (i in 1:iter) {                            
      #y <- unlist(mapply(generate, n=cluster.sizes, mu=cond_means, ...))
      y_matrix[i, ] <- generate()#y
  }

  result <- structure(list("x" = x, "z" = z, "y" = y_matrix, "family" = family, "cluster" = cluster_assignments), class = "heritData")
  return(result)
}

print.heritData <- function(dataObj) {
cat("Family:", dataObj$family, "\n")
  
  nx <- (ncol(dataObj$x)) - 1
  datamat <- rbind(t(dataObj$x), dataObj$y)
  colnames(datamat) <- dataObj$cluster
  rownames(datamat) <- c(paste0('x', 0:nx), paste0('gene', 1:nrow(dataObj$y)))
  
  print(datamat)
}

# Example
x <- c(23,24,25,26,29,30,18,22,21,16,18,16,21)
beta <- c(7,9)
sigma.u <- matrix(c(2,1,1,2),2)
cluster.sizes <- c(3,2,5,3)
theta = 2
phi = 2
power = 1.6
nu = 3
lambda2 = 0.7
iter =10
data1 <- herit.glmmdata(iter, x, beta, sigma.u, cluster.sizes, "nbinom2", theta=theta)
data2 <- herit.glmmdata(iter, x, beta, sigma.u, cluster.sizes, "tweedie", phi=phi, power=power)
data3 <- herit.glmmdata(iter, x, beta, sigma.u, cluster.sizes, "compois", nu=nu)
#data4 <- herit.glmmdata(iter, x, beta, sigma.u, cluster.sizes, "genpois", lambda2=lambda2)
