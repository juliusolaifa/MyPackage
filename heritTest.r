herit.test <- function(object, ...) {
  UseMethod("herit.test", object)
}

herit.test.heritScore <- function(scoreObj, new_formula) {
  modObj <- scoreObj$modObj
  return(herit.test.heritMod(modObj,new_formula))
}

herit.test.heritMod <- function(modObj, new_formula) {
  modObj_extract <- modObj$modObj
  L_H1 <- modObj$logLikelihood
  total_par <- modObj$npar
  total_varcomp <- modObj$npar - modObj$nfixef
  npar_bounded <- floor((sqrt(1 + 8 * total_varcomp) - 1) / 2) # S = n(n+1)/2, where S =n(diag)
    
  reduced_modObj <- update(modObj_extract, formula = new_formula)
  L_H0 <- as.numeric(logLik(reduced_modObj))
  total_reduced_npar <- attr(logLik(reduced_modObj), "df")
  total_reduced_nfixef <- length(c(fixef(reduced_modObj)$cond, sigma(reduced_modObj)))
  total_reduced_varcomp <- total_reduced_npar - total_reduced_nfixef
  npar_reduced_bounded <- floor((sqrt(1 + 8 * total_reduced_varcomp) - 1) / 2)

  nvarcomp_on_boundary <- npar_bounded - npar_reduced_bounded
  likStat <- -2*(L_H0 - L_H1)

  ###
  critical_value <- compute_chisq_critical()
  p_value <- 
  # Return the results
  return(list(test_results = "Results of your test"))
}

compute_chisq_critical <- function() {
}
