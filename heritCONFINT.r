confint.heritScore <- function(object, parm, level = 0.95, ...) {
  # Extract the heritability estimate from the object
  heritability <- object$heritability
  
  # Calculate the confidence interval (this is just an example; you would replace this with the appropriate calculation based on your model)
  ci_low <- heritability - 1.96 * sqrt(heritability * (1 - heritability) / n)
  ci_high <- heritability + 1.96 * sqrt(heritability * (1 - heritability) / n)
  
  # Create a data frame to store the results
  result <- data.frame(
    Estimate = heritability,
    `2.5 %` = ci_low,
    `97.5 %` = ci_high
  )
  
  return(result)
}

confint.heritMod <- function(object, parm, level = 0.95, ...) {
  # Extract the model parameters from the object
  beta <- object$beta
  Sigma <- object$Sigma
  phi <- object$phi
  
  # Calculate the confidence intervals for the parameters (this is just an example; you would replace this with the appropriate calculations based on your model)
  beta_ci <- confint(object$modObj, parm = "beta", level = level)
  Sigma_ci <- matrix(c(confint(object$modObj, parm = "Sigma", level = level)), ncol = 2)
  phi_ci <- c(phi - 1.96 * sqrt(phi), phi + 1.96 * sqrt(phi))
  
  # Create a list to store the results
  result <- list(
    beta = beta_ci,
    Sigma = Sigma_ci,
    phi = phi_ci
  )
  
  return(result)
}
