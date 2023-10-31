herit.test <- function(object, ...) {
  UseMethod("herit.test", object)
}

herit.test.heritScore <- function(scoreObj, new_formula) {
  modObj <- scoreObj$modObj
  return(herit.test.heritMod(modObj,new_formula))
}

herit.test.heritMod <- function(modObj, new_formula) {
  
  full_modObj <- modObj$modObj
  full_logLik <- full_modObj$logLikelihood
  full_nfixef <- modObj$nfixef
  full_nvarcomp_bounded <- modObj$nvarcomp_bounded
    
  reduced_modObj <- update(full_modObj, formula = new_formula)
  reduced_logLik <- as.numeric(logLik(reduced_modObj))
  reduced_nfixef <- reduced_modObj$nfixef
  reduced_nvarcomp_bounded <- reduced_modObj$nvarcomp_bounded

  likStat <- -2*(L_H0 - L_H1)

  ###
  pvalue_value <- compute_chisq_pvalue(full_nfixef, reduced_nfixef, full_nvarcomp_bounded, reduced_nvarcomp_bounded)
  # Return the results
  return(list(test_results = "Results of your test"))
}

compute_chisq_pvalue <- function(full_nfixef, reduced_nfixef, full_nvarcomp_bounded, reduced_nvarcomp_bounded) {
    if(full_nfixef == reduced_nfixef) {
      
    }else if(full_nfixef = reduced_nfixef + 1) {
      
    }else{
      stop("This case is not implemented)
    }
}
