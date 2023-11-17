# obtain heritability score
herit.score <- function(data, formula, X, Z=NULL, family = nbinom2, REML = TRUE) {
  fitmod <- herit.fit(data, formula, family = family, REML = REML)
  herit <- herit.vpc(fitmod, X, Z)
  return(structure(list("heritability"=herit, "nfixef"=fitmod$nfixef, "nvarcomp_bounded"=fitmod$nvarcomp_bounded
                        "logLikelihood"=logLike,"modObj"=fitmod$modObj), class="heritScore"))
}

# fit a glmm
herit.fit <- function(fixed, random, data, family) {
    modObj <- GLMMadaptive::mixed_model(fixed, random, data, family)
    # groupVar <- modObj$modelInfo$grpVar
    # link <- modObj$modelInfo$family$link
    # beta <- fixef(modObj)$cond
    # std <- attr(VarCorr(modObj)[[c("cond", groupVar)]], "stddev")
    # corr <- attr(VarCorr(modObj)[[c("cond", groupVar)]], "correlation")
    # phi <- sigma(modObj)
    # Sigma <- corr*(std %*% std)
    # nfixef <- length(c(beta,phi))
    # nvarcomp_bounded <- length(std)
    # logLike <- as.numeric(logLik(modObj))
    # result <- list("beta"=beta, "Sigma"=Sigma, "phi"=phi, "family"=family, 
    #                "link"=link, "nfixef"=nfixef, "nvarcomp_bounded"=nvarcomp_bounded, "logLikelihood"=logLike, "modObj"=modObj)
    # class(result) <- "heritMod"
    # return(result)
}

# Define the generic function
herit.vpc <- function(X, Z=NULL, beta, Sigma, phi, family, link) {
  UseMethod("herit.vpc")
}

# Define the default method (direct parameter input)
herit.vpc.default <- function(X,Z=NULL,beta,Sigma,phi,family,link) {
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
    return(vpc_compute(mu, sigm, phi, family, link))
}

# Define the method for heritMod input
herit.vpc.heritMod <- function(modelObj, X, Z=NULL) {
  # Extract parameters from the fitted model object
  beta <- modelObj$beta
  Sigma <- modelObj$Sigma
  phi <- modelObj$phi
  family <- modelObj$family
  link <- modelObj$link
  
  # Call the default method with extracted parameters
  herit.vpc.default(X, Z, beta, Sigma, phi, family, link)
}

# this methods computes the actual vpc for a specified family using the right variance mean function
vpc_compute <- function(mu, sigm, phi, family, link, p=NULL) {

    if(link == "log") {
        inv.mu <- exp(mu + sigm/2)
        inv.var <- (exp(sigm) - 1)*exp(2*mu + sigm)
    }
    switch(family,
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
