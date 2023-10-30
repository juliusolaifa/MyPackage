# Define the generic function
herit.vpc <- function(X, Z=NULL, beta, Sigma, phi, family, g_inv.dist) {
  UseMethod("herit.vpc")
}

# Define the default method (direct parameter input)
herit.vpc.default <- function(X,Z=NULL,beta,Sigma,phi,family,g_inv.dist) {
    if(is.null(Z)) {
        Z = X
    }

    # Check if X and beta have compatible dimensions
    if (length(X) != length(beta)) {
        stop("Dimensions of X and beta are not compatible for multiplication.")
    }

    # Check if X and beta have compatible dimensions
    if (length(Z) != ncol(Sigma)) {
        stop("Dimensions of Z and Sigma are not compatible for matrix multiplication.")
    }
    
    mu <- X%*%beta
    sigm <- Z%*%Sigma%*%Z
    return(vpc_compute(mu, sigm, phi, family, g_inv.dist))
}

# Define the method for fitted_model input
herit.vpc.fitted_model <- function(modelObj, X, Z=NULL) {
  # Extract parameters from the fitted model object
  beta <- modelObj$coefficients$beta
  Sigma <- modelObj$coefficients$Sigma
  phi <- modelObj$coefficients$phi
  family <- modelObj$family
  g_inv.dist <- modelObj$g_inv.dist
  
  # Call the default method with extracted parameters
  herit.vpc.default(X, Z, beta, Sigma, phi, family, g_inv.dist)
}

# this methods computes the actual vpc for a specified family using the right variance mean function
vpc_compute <- function(mu, sigm, phi, g_inv.dist, family, p=NULL) {

    if(g_inv.dist == "log") {
        inv.mu <- exp(mu + sigm/2)
        inv.var <- (exp(sigm) - 1)*exp(2*mu + sigm)
    }
    switch(family,
           nbimom1={
               return(inv.var / (inv.var + inv.mu + phi*inv.mu))
           }
           nbinom2={
               inv.mu.p <- exp(2*mu + 4*sigm/2)
               return(inv.var / (inv.var + inv.mu + inv.mu.p/phi))
           },
           tweedie={
               if(is.null(p) || p <= 1 || p >= 2) {
                   stop("p must be an integer between 1 and 2")
               }
               inv.mu.p <- exp(p*mu + p^2*sigm/2)
               return(inv.var / (inv.var + phi*inv.mu.p))
           }
           compois={
               return(inv.var / (inv.var + phi*inv.mu))
           }
           genpois={
               return(inv.var / (inv.var + phi*inv.mu)) #phi = exp(eta)
           }
           {
               return(paste("VPC not implemented for family: ", family))
           }
    )
}

fitted_model_obj <- structure(list(), class = "fitted_model")
params_obj <- structure(list(), class = "params")
