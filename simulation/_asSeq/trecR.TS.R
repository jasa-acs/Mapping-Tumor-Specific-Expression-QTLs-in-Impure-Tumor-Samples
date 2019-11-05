
#-------------------------------------------------------------
# log likelihood of Negative-Binomial
#-------------------------------------------------------------

loglikNB <- function(y, muTReC, phi){
  
  logL = 0.0
  
  vphi = 1/phi
  lik0 = vphi*log(vphi) - lgamma(vphi)
  
  for(i in 1:length(y)){
    yi   = y[i]
    mui  = muTReC[i]
    
    if(yi==0){
      logL = logL + vphi*log(vphi) - vphi*log(vphi + mui)
    }else{
      logL = logL + lgamma(yi + vphi) - lgamma(yi + 1.0) + yi*log(mui)
      logL = logL - (vphi+yi)*log(vphi+mui) + lik0
    }
  }
  
  logL
}

# -----------------------------------------------------------------
# likelihood of TReC data with respect to log(kappa), log(eta),
# and log(gamma)
#
# here mu0 = exp(Xb), where X are all the covariates
# (including intercept) other than kappa, eta, and gamma
# -----------------------------------------------------------------

neg.logLTReC.TS <- function(para, y, z1, nu0, phi, rhos, H0, gH=FALSE){
  
  kappa = exp(para[1])
  
  if(H0==""){
    eta   = exp(para[2])
    gamma = exp(para[3])
    
  }else if(H0=="eta"){
    eta   = 1.0
    gamma = exp(para[2])
    
  }else if(H0=="gamma"){
    eta   = exp(para[2])
    gamma = 1.0
    
  }else{
    stop("invalid value for H0!\n")
  }
  
  delta  = (1-rhos) + rhos*kappa
  cs     = (1-rhos)/delta
  zeta   = cs*eta + (1-cs)*gamma

  ww1    = which(z1==1)
  ww2    = which(z1==2)

  nu1      = nu0 * delta
  nu1[ww1] = nu1[ww1] * (1 + zeta[ww1])/2
  nu1[ww2] = nu1[ww2] * zeta[ww2]

  logL = -loglikNB(y, nu1, phi)
  
  if(gH){
    
    gh1 = neg.gradTReC.TS(para, y, z1, nu0, phi, rhos, H0, Hessian=TRUE)
    
    attr(logL, "gradient") = gh1$gr
    attr(logL, "hessian")  = gh1$hs
  }
  
  logL
}

# -----------------------------------------------------------------
# gradient and Hessian of TReC log likelihood with respect to
# log(kappa), log(eta), and log(gamma)
# -----------------------------------------------------------------

neg.gradTReC.TS <- function(para, y, z1, nu0, phi, rhos, H0, Hessian=FALSE){
  
  kappa = exp(para[1])
  
  if(H0==""){
    eta   = exp(para[2])
    gamma = exp(para[3])
    
  }else if(H0=="eta"){
    eta   = 1.0
    gamma = exp(para[2])
    
  }else if(H0=="gamma"){
    eta   = exp(para[2])
    gamma = 1.0
    
  }else{
    stop("invalid value for H0!\n")
  }
  
  delta  = (1-rhos) + rhos*kappa
  cs     = (1-rhos)/delta
  zeta   = cs*eta + (1-cs)*gamma
  
  ww0    = which(z1==0)
  ww1    = which(z1==1)
  ww2    = which(z1==2)
  
  nu1      = nu0*delta
  nu1[ww1] = nu1[ww1]*(1 + zeta[ww1])/2
  nu1[ww2] = nu1[ww2]*zeta[ww2]

  dl.dmu   = y/nu1 - (1 + phi*y)/(1 + phi*nu1)
  
  # ---------------------------------------------
  # d_logL/d_log_kappa
  # ---------------------------------------------
  
  dmu.dkappa      = rep(NA, length(nu0))
  dmu.dkappa[ww0] = nu0[ww0] * rhos[ww0]
  dmu.dkappa[ww1] = 0.5*nu0[ww1] * rhos[ww1] * (1 + gamma)
  dmu.dkappa[ww2] = nu0[ww2] * rhos[ww2] * gamma
  
  dl.dkappa  = kappa*sum(dl.dmu * dmu.dkappa)
  
  # ---------------------------------------------
  # d_logL/d_log_eta
  # ---------------------------------------------
  
  if(H0=="eta"){
    dl.deta = 0.0
  }else{
    dmu.deta      = rep(0, length(nu0))
    dmu.deta[ww1] = 0.5*nu0[ww1]*(1 - rhos[ww1])
    dmu.deta[ww2] = nu0[ww2]*(1 - rhos[ww2])
    
    dl.deta = eta*sum(dl.dmu * dmu.deta)
  }
  
  # ---------------------------------------------
  # d_logL/d_log_gamma
  # ---------------------------------------------

  if(H0=="gamma"){
    dl.dgamma = 0.0
  }else{
    dmu.dgamma      = rep(0, length(nu0))
    dmu.dgamma[ww1] = 0.5*nu0[ww1]*rhos[ww1]*kappa
    dmu.dgamma[ww2] = nu0[ww2]*rhos[ww2]*kappa
    
    dl.dgamma = gamma*sum(dl.dmu * dmu.dgamma)
  }
  
  if(H0==""){
    grd = c(dl.dkappa, dl.deta, dl.dgamma)
  }else if(H0=="eta"){
    grd = c(dl.dkappa, dl.dgamma)
  }else if(H0=="gamma"){
    grd = c(dl.dkappa, dl.deta)
  }

  #-----------------------------------------------------------
  # calculate Hessian matrix of TReC likelihood with respect
  # to log(kappa), log(eta), and log(gamma)
  #-----------------------------------------------------------

  if(Hessian){
    d2l.dmu2  = -y/(nu1*nu1) + phi*(1 + phi*y)/((1 + phi*nu1)^2)
    
    sum1 = sum(d2l.dmu2 * dmu.dkappa * dmu.dkappa)
    d2l.dkappa2 = dl.dkappa + kappa * kappa * sum1
    
    if(H0!="eta"){
      sum1 = sum(d2l.dmu2 * dmu.deta * dmu.deta)
      d2l.deta2 = dl.deta + eta * eta * sum1
      
      d2l.dkappa.deta = eta * kappa * sum(d2l.dmu2 * dmu.deta * dmu.dkappa)
    }

    if(H0!="gamma"){
      sum1 = sum(d2l.dmu2 * dmu.dgamma * dmu.dgamma)
      d2l.dgamma2 = dl.dgamma + gamma * gamma * sum1
      
      d2mu.dkappa.dgamma = rep(0, length(nu0))
      d2mu.dkappa.dgamma[ww1] = 0.5*nu0[ww1]*rhos[ww1]
      d2mu.dkappa.dgamma[ww2] = nu0[ww2]*rhos[ww2]
      
      sum1 = sum(d2l.dmu2 * dmu.dgamma * dmu.dkappa)
      sum2 = sum(dl.dmu * d2mu.dkappa.dgamma)
      d2l.dkappa.dgamma = gamma * kappa * (sum1 + sum2)
    }
    
    if(H0==""){
      d2l.deta.dgamma = eta * gamma * sum(d2l.dmu2 * dmu.dgamma * dmu.deta)
      
      hes = matrix(NA, nrow=3, ncol=3)
      diag(hes) = c(d2l.dkappa2, d2l.deta2, d2l.dgamma2)
      hes[1,2] = hes[2,1] = d2l.dkappa.deta
      hes[1,3] = hes[3,1] = d2l.dkappa.dgamma
      hes[2,3] = hes[3,2] = d2l.deta.dgamma
      
    }else if(H0 == "eta"){
      hes = matrix(NA, nrow=2, ncol=2)
      diag(hes) = c(d2l.dkappa2, d2l.dgamma2)
      hes[1,2]  = hes[2,1] = d2l.dkappa.dgamma
      
    }else if(H0 == "gamma"){
      hes = matrix(NA, nrow=2, ncol=2)
      diag(hes) = c(d2l.dkappa2, d2l.deta2)
      hes[1,2] = hes[2,1] = d2l.dkappa.deta
    }
  }
  
  if(Hessian){
    res = list(gr=-grd, hs=-hes)
  }else{
    res = -grd
  }
  
  res
}

# -----------------------------------------------------------------
# Hessian of TReC log likelihood with respect to
# log(kappa), log(eta), and log(gamma)
# -----------------------------------------------------------------

neg.HessianTReC.TS <- function(para, y, z1, nu0, phi, rhos, H0){
  
  kappa = exp(para[1])
  
  if(H0==""){
    eta   = exp(para[2])
    gamma = exp(para[3])
    
  }else if(H0=="eta"){
    eta   = 1.0
    gamma = exp(para[2])
    
  }else if(H0=="gamma"){
    eta   = exp(para[2])
    gamma = 1.0
    
  }else{
    stop("invalid value for H0!\n")
  }
  
  delta  = (1-rhos) + rhos*kappa
  cs     = (1-rhos)/delta
  zeta   = cs*eta + (1-cs)*gamma
  
  ww0    = which(z1==0)
  ww1    = which(z1==1)
  ww2    = which(z1==2)
  
  nu1      = nu0*delta
  nu1[ww1] = nu1[ww1]*(1 + zeta[ww1])/2
  nu1[ww2] = nu1[ww2]*zeta[ww2]
  
  dl.dmu   = y/nu1 - (1 + phi*y)/(1 + phi*nu1)
  
  # ---------------------------------------------
  # d_logL/d_log_kappa
  # ---------------------------------------------
  
  dmu.dkappa      = rep(NA, length(nu0))
  dmu.dkappa[ww0] = nu0[ww0] * rhos[ww0]
  dmu.dkappa[ww1] = 0.5*nu0[ww1] * rhos[ww1] * (1 + gamma)
  dmu.dkappa[ww2] = nu0[ww2] * rhos[ww2] * gamma
  
  dl.dkappa  = kappa*sum(dl.dmu * dmu.dkappa)
  
  # ---------------------------------------------
  # d_logL/d_log_eta
  # ---------------------------------------------
  
  if(H0=="eta"){
    dl.deta = 0.0
  }else{
    dmu.deta      = rep(0, length(nu0))
    dmu.deta[ww1] = 0.5*nu0[ww1]*(1 - rhos[ww1])
    dmu.deta[ww2] = nu0[ww2]*(1 - rhos[ww2])
    
    dl.deta = eta*sum(dl.dmu * dmu.deta)
  }
  
  # ---------------------------------------------
  # d_logL/d_log_gamma
  # ---------------------------------------------
  
  if(H0=="gamma"){
    dl.dgamma = 0.0
  }else{
    dmu.dgamma      = rep(0, length(nu0))
    dmu.dgamma[ww1] = 0.5*nu0[ww1]*rhos[ww1]*kappa
    dmu.dgamma[ww2] = nu0[ww2]*rhos[ww2]*kappa
    
    dl.dgamma = gamma*sum(dl.dmu * dmu.dgamma)
  }
  
  #-----------------------------------------------------------
  # calculate Hessian matrix of TReC likelihood with respect
  # to log(kappa), log(eta), and log(gamma)
  #-----------------------------------------------------------
  
  d2l.dmu2  = -y/(nu1*nu1) + phi*(1 + phi*y)/((1 + phi*nu1)^2)
  
  sum1 = sum(d2l.dmu2 * dmu.dkappa * dmu.dkappa)
  d2l.dkappa2 = dl.dkappa + kappa * kappa * sum1
  
  if(H0!="eta"){
    sum1 = sum(d2l.dmu2 * dmu.deta * dmu.deta)
    d2l.deta2 = dl.deta + eta * eta * sum1
    
    d2l.dkappa.deta = eta * kappa * sum(d2l.dmu2 * dmu.deta * dmu.dkappa)
  }
  
  if(H0!="gamma"){
    sum1 = sum(d2l.dmu2 * dmu.dgamma * dmu.dgamma)
    d2l.dgamma2 = dl.dgamma + gamma * gamma * sum1
    
    d2mu.dkappa.dgamma = rep(0, length(nu0))
    d2mu.dkappa.dgamma[ww1] = 0.5*nu0[ww1]*rhos[ww1]
    d2mu.dkappa.dgamma[ww2] = nu0[ww2]*rhos[ww2]
    
    sum1 = sum(d2l.dmu2 * dmu.dgamma * dmu.dkappa)
    sum2 = sum(dl.dmu * d2mu.dkappa.dgamma)
    d2l.dkappa.dgamma = gamma * kappa * (sum1 + sum2)
  }
  
  if(H0==""){
    d2l.deta.dgamma = eta * gamma * sum(d2l.dmu2 * dmu.dgamma * dmu.deta)
    
    hes = matrix(NA, nrow=3, ncol=3)
    diag(hes) = c(d2l.dkappa2, d2l.deta2, d2l.dgamma2)
    hes[1,2] = hes[2,1] = d2l.dkappa.deta
    hes[1,3] = hes[3,1] = d2l.dkappa.dgamma
    hes[2,3] = hes[3,2] = d2l.deta.dgamma
    
  }else if(H0 == "eta"){
    hes = matrix(NA, nrow=2, ncol=2)
    diag(hes) = c(d2l.dkappa2, d2l.dgamma2)
    hes[1,2]  = hes[2,1] = d2l.dkappa.dgamma
    
  }else if(H0 == "gamma"){
    hes = matrix(NA, nrow=2, ncol=2)
    diag(hes) = c(d2l.dkappa2, d2l.deta2)
    hes[1,2] = hes[2,1] = d2l.dkappa.deta
  }
  
  res = -hes
  
  res
}

#-------------------------------------------------------------
# test plot, plot gradient and likelihood across the range
# of one of the parameters
#-------------------------------------------------------------

plotItTReC <- function(y, z1, nu0, phi, rhos, kappa, eta, gamma){
  
  ## tmp are the values of that parameter to be checked will take
  nk   = 100
  tmp  = seq(exp(-1.5), exp(1.5), length.out=nk)
  
  ## check for kappa, gamma, and eta
  logL = grad = list()
  
  for(kk in 1:3){
    
    if(kk == 1){
      kappas = tmp
      etas   = rep(eta, nk)
      gammas = rep(gamma, nk)
    }
    
    if(kk == 2){
      kappas = rep(kappa, nk)
      etas   = tmp
      gammas = rep(gamma,   nk)
    }
    
    if(kk == 3){
      kappas = rep(kappa, nk)
      etas   = rep(eta,   nk)
      gammas = tmp
    }
    
    ll = rep(NA, nk)
    gs = matrix(NA, nrow=nk, ncol=3)
    
    for(i in 1:nk){
      para   = log(c(kappas[i], etas[i], gammas[i]))
      ll[i]  = -neg.logLTReC.TS(para, y, z1, nu0, phi, rhos, H0="")
      gs[i,] = -neg.gradTReC.TS(para, y, z1, nu0, phi, rhos, H0="")
    }
    
    logL[[kk]] = ll
    grad[[kk]] = gs
  }
  
  nms = c("kappa", "eta", "gamma")
  vls = list(kappa, eta, gamma)
  names(vls) = nms
  
  par(mfrow=c(2,3), mar=c(5,4,2,1), bty="n")
  
  for(kk in 1:3){
    plot(tmp, grad[[kk]][,kk], type="l", xlab=nms[kk], ylab="gradient")
    abline(h=0, lty=2)
    abline(v=1, col="grey")
    abline(v=vls[[nms[kk]]], col="red")
  }
  
  for(kk in 1:3){
    plot(tmp, logL[[kk]], type="l", xlab=nms[kk], ylab="log-likelihood")
    abline(v=1, col="grey")
    abline(v=vls[[nms[kk]]], col="red")
  }
  
}

#-------------------------------------------------------------
# TReC model
#-------------------------------------------------------------

trecR.TS <- function(y, Xs, z1, rhos, H0, para0=NULL, nIter=500,
                     traceIt=FALSE, method=c("BFGS", "Newton", "nlminb"),
                     yfit=FALSE, convergence=1e-4){
  
  if(! all(unique(z1) %in% c(0,1,2))){
    stop("z1 has invlaid values\n")
  }
  
  ww1  = which(z1==1)
  ww2  = which(z1==2)

  gc1  = glm.control(epsilon=1e-5, maxit=25, trace = FALSE)
  
  #----------------------------------------------------------
  # initial model fitting
  #----------------------------------------------------------
  
  if(! is.null(para0)){
    kappa = para0[1]
    
    if(H0==""){
      eta   = para0[2]
      gamma = para0[3]
      
    }else if(H0=="eta"){
      eta   = 1.0
      gamma = para0[2]
      
    }else if(H0=="gamma"){
      eta   = para0[2]
      gamma = 1.0
      
    }else{
      stop("invalid value for H0!\n")
    }
    
    #-----------------------------------------------------------
    # estimate b and phi
    #-----------------------------------------------------------
    
    off1 = log(1-rhos + rhos*kappa)
    
    cs   = (1-rhos)/(1-rhos + rhos*kappa)
    zeta = cs * eta + (1 - cs)*gamma
    off1[ww1] = off1[ww1] + log((zeta[ww1] + 1)/2)
    off1[ww2] = off1[ww2] + log(zeta[ww2])
    
    g0  = glm.nb(y ~ -1 + Xs + offset(off1),  control=gc1)

    nu0   = exp(Xs %*% g0$coef)
    phi0  = 1/g0$theta

  }else{
    g0    = glm.nb(y ~ -1 + Xs, control=gc1)
    
    nu0   = exp(Xs %*% g0$coef)
    phi0  = 1/g0$theta
    
    kappa = 1.0
    eta   = 1.0
    gamma = 1.0
    
    if(H0==""){
      para0 = c(kappa, eta, gamma)
      names(para0) = c("kappa", "eta", "gamma")
      
    }else if(H0=="eta"){
      para0 = c(kappa, gamma)
      names(para0) = c("kappa", "gamma")
      
    }else if(H0=="gamma"){
      para0 = c(kappa, eta)
      names(para0) = c("kappa", "eta")
    }
  }
  
  logLik  = NULL
  logLik0 = -Inf

  if(traceIt)
  message(sprintf("Initial: kappa=%.2e, eta=%.2e, gamma=%.2e, phi0=%.2e\n",
                  kappa, eta, gamma, phi0))

  for(g in 1:nIter){
    
    parDiff = 0.0

    #--------------------------------------------------------
    # estimate parametrs first
    #--------------------------------------------------------

    if(method=="BFGS"){
      
      o1 = optim(log(para0), fn=neg.logLTReC.TS, gr=neg.gradTReC.TS, y=y,
                z1=z1, nu0=nu0, phi=phi0, rhos=rhos, H0=H0,
                method="BFGS", control=list(maxit=10, trace=0))
      
      # o1$convergence = 0 means success
      # and o1$convergence = 1 means reaching iteration limit
      if(o1$convergence != 0 && o1$convergence != 1){
        warning("g =", g, " fail BFGS\n")
      }
      
      para1   = exp(o1$par)
      logLik1 = -o1$value

    }else if(method=="Newton"){
      
      n1 = nlm(f=neg.logLTReC.TS, p=log(para0), y=y, z1=z1, nu0=nu0,
                phi=phi0, rhos=rhos, H0=H0, gH=TRUE, iterlim=10)
      # o1$convergence = 0 means success
      # and o1$convergence = 1 means reaching iteration limit
      if(n1$code == 3 || n1$code == 5){
        warning("g =", g, " nlm may not converge.\n")
      }

      para1   = exp(n1$estimate)
      names(para1) = names(para0)
      logLik1 = -n1$minimum
      
    }else if(method=="nlminb"){
      
      nl1 = nlminb(log(para0), neg.logLTReC.TS, gradient = neg.gradTReC.TS,
                  hessian=neg.HessianTReC.TS, y=y, z1=z1, nu0=nu0, phi=phi0,
                  rhos=rhos, H0=H0, control = list(iter.max=10))
                  
      para1   = exp(nl1$par)
      logLik1 = -nl1$objective

    }else{
      stop("invalid value for method\n.")
      
    }
    
    if(traceIt){
      message(sprintf("\ng=%d", g))
      message(sprintf("  Estimate para: logLik_TReC=(%.3e, %.3e)",
              logLik0, logLik1))
    }
    
    if(g > 1){
      if((logLik1 - logLik0)/abs(logLik0) < -1e-4){
        stop("liklihood decreases for para\n")
      }
    }
    
    if(parDiff < max(abs(para1 - para0))){
      parDiff = max(abs(para1 - para0))
    }
    
    para0    = para1
    logLik0  = logLik1
    
    kappa = para0[1]
    
    if(H0==""){
      eta   = para0[2]
      gamma = para0[3]
      
    }else if(H0=="eta"){
      eta   = 1.0
      gamma = para0[2]
      
    }else if(H0=="gamma"){
      eta   = para0[2]
      gamma = 1.0
      
    }else{
      stop("invalid value for H0!\n")
    }
    
    #-----------------------------------------------------------
    # estimate b and phi
    #-----------------------------------------------------------
    
    off1 = log(1-rhos + rhos*kappa)
    
    cs   = (1-rhos)/(1-rhos + rhos*kappa)
    zeta = cs * eta + (1 - cs)*gamma
    off1[ww1] = off1[ww1] + log((zeta[ww1] + 1)/2)
    off1[ww2] = off1[ww2] + log(zeta[ww2])
    
    g1      = glm.nb(y ~ -1 + Xs + offset(off1),  control=gc1)
    phi1    = 1/g1$theta
    nu1off  = exp(Xs %*% g1$coef + off1)
    logLik1 = loglikNB(y, nu1off, phi1)
    
    if((logLik1 - logLik0)/abs(logLik0) < -1e-5){
      stop("liklihood decreases for estimating b and phi\n")
    }
    
    if(parDiff < abs(phi1 - phi0)){
      parDiff = abs(phi1 - phi0)
    }
    
    if(traceIt){
      message(sprintf("Estimate b & phi: phi=%.2e, logLikTReC=(%.3e, %.3e)",
            phi1, logLik0, logLik1))
    }
    
    phi0 = phi1
    nu0  = exp(Xs %*% g1$coef)
    logLik0 = logLik1

    #-----------------------------------------------------------
    # check convergence
    #-----------------------------------------------------------

    logLik = c(logLik, logLik1)
    
    if(traceIt){
      str = "g=%d, parDiff=%.2e, kappa=%.2f, eta=%.2f, gamma=%.2f, phi0=%.2e"
      message(sprintf(str, g, parDiff, kappa, eta, gamma, phi0))
    }
    
    if(parDiff < convergence) break
    
  }
    
  l1 = list(para=para0, phi=phi0, bs=g1$coef, logLik=logLik)

  if(yfit){
    l1$fitted = g1$fitted
  }
  
  l1
}

#-------------------------------------------------------------
# TReC model
# z1 is genotype data coded as 0, 1, and 2
#-------------------------------------------------------------

trecR.TS.test <- function(y, Xs, z1, rhos, nIter=500,
                          method=c("BFGS", "Newton", "nlminb"),
                          traceIt=FALSE)
{
  
  t1      = trecR.TS(y, Xs, z1, rhos, nIter=nIter, H0="",
                     method=method, traceIt=traceIt)
  
  t0eta   = trecR.TS(y, Xs, z1, rhos, nIter=nIter, H0="eta",
                     para0=t1$para[c(1,3)], method=method,
                     traceIt=traceIt)
  
  t0gamma = trecR.TS(y, Xs, z1, rhos, nIter=nIter, H0="gamma",
                     para0=t1$para[c(1,2)], method=method,
                     traceIt=traceIt)
  
  l0eta   = t0eta$logLik[length(t0eta$logLik)]
  l0gamma = t0gamma$logLik[length(t0gamma$logLik)]
  l1      = t1$logLik[length(t1$logLik)]
  
  lrt.eta   = 2*(l1 - l0eta)
  lrt.gamma = 2*(l1 - l0gamma)
  
  pv.eta    = pchisq(lrt.eta, df=1, lower.tail=FALSE)
  pv.gamma  = pchisq(lrt.gamma, df=1, lower.tail=FALSE)
  
  pvs = data.frame(LRS=c(lrt.eta, lrt.gamma), pval=c(pv.eta, pv.gamma))
  rownames(pvs) = c("eta", "gamma")
  
  list(para=t1$para, phi=t1$phi, bs=t1$bs, pvalue=pvs)
}

