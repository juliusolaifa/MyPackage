vcov.VarCorr.merMod <- function(object,fit,...) {
    
    if (isREML(fit)) {
        warning("refitting model with ML")
        fit <- refitML(fit)
    }
    if (!require("numDeriv")) stop("numDeriv package required")
    useSc <- attr(object,"useSc")
    print(useSc)
    dd <- lme4:::devfun2(fit,useSc=useSc,signames=FALSE)
    vdd <- as.data.frame(object,order="lower.tri")
    pars <- vdd[,"sdcor"]
    npar0 <- length(pars)
    if (isGLMM(fit)) {
        pars <- c(pars,lme4::fixef(fit))
    }
    hh1 <- hessian(dd,pars)
    vv2 <- 2*solve(hh1)
    if (isGLMM(fit)) {
        vv2 <- vv2[1:npar0,1:npar0,drop=FALSE]
    }
    nms <- apply(vdd[,1:3],1,
                 function(x) paste(na.omit(x),collapse="."))
    dimnames(vv2) <- list(nms,nms)
    return(vv2)
}