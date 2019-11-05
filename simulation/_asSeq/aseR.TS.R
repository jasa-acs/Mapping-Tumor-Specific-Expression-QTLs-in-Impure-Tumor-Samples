
#-------------------------------------------------------------
# log likelihood of Beta-Binomial
#-------------------------------------------------------------

loglikBB <- function(nB, nTotal, pis, psi){
  
  if(any(nTotal < 1)){
    stop("nTotal must be positive\n")
  }
  
  logL = 0.0
  
  for(i in 1:length(nTotal)){
    
    ni   = nTotal[i]
    nBi  = nB[i]
    pi1  = pis[i]
    
    logLi = lchoose(ni, nBi)
    
    if(nBi > 0){
      k     = 0:(nBi-1)
      logLi = logLi + sum(log(pi1 + k*psi))
    }
    
    if(nBi < ni){
      k     = 0:(ni-nBi-1)
      logLi = logLi + sum(log(1 - pi1 + k*psi))
    }
    
    if(ni > 1){
      k     = 1:(ni-1)
      logLi = logLi - sum(log(1 + k*psi))
    }

    logL  = logL + logLi
  }
  
  logL
}

# -----------------------------------------------------------------
# log likelihood of ASE data with respect to log(kappa), log(eta),
# and log(gamma)
# -----------------------------------------------------------------

neg.logLASE.TS <- function(para, nB, nTotal, psi, rhos, H0, gH=FALSE){
  
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
  pis   = zeta/(zeta + 1)

  logL = -loglikBB(nB, nTotal, pis, psi)
  
  if(gH){
    
    gh1 = neg.gradASE.TS(para, nB, nTotal, psi, rhos, H0, Hessian=TRUE)
    
    attr(logL, "gradient") = gh1$gr
    attr(logL, "hessian")  = gh1$hs
  }
  
  logL
}

# --------------------------------------------------------------    
# Gradient of logLik of beta-binomial with respect to psi
# --------------------------------------------------------------

gradASE.psi <- function(psi, nB, nTotal, pis){
  
  grPsi = 0
  
  for(i in 1:length(nB)){
    
    ni   = nTotal[i]
    nBi  = nB[i]
    pi   = pis[i]
    
    if(nBi > 0){
      k     = seq(0, nBi-1, by=1)
      grPsi = grPsi + sum(k/(pi + k*psi))
    }
    
    if(nBi < ni){
      k     = seq(0, ni-nBi-1, by=1)
      grPsi = grPsi + sum(k/(1 - pi + k*psi))
    }
    
    if(ni > 1){
      k     = 1:(ni-1)
      grPsi = grPsi - sum(k/(1 + k*psi))
    }
  }
  
  grPsi
}

# --------------------------------------------------------------
# Gradient of logLik of beta-binomial with respect to
# log(kappa), log(eta), and log(gamma)
#
# here we assume all the samples are those with heterozygous
# genotype at the candidate eQTL
# --------------------------------------------------------------

neg.gradASE.TS <- function(para, nB, nTotal, psi, rhos, H0, Hessian=FALSE){
    
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
  pis    = zeta/(1 + zeta)

  dpi.dzeta = 1/((1 + zeta)^2)
  
  # ---------------------------------------------
  # dl.dpi
  # ---------------------------------------------

  dl.dpi = rep(0, length(nB))

  for(i in 1:length(nB)){
    
    grPi = 0

    ni   = nTotal[i]
    nBi  = nB[i]
    pi   = pis[i]
    
    if(nBi > 0){
      k    = seq(0, nBi-1, by=1)
      grPi = grPi + sum(1/(pi + k*psi))
    }
    
    if(nBi < ni){
      k    = seq(0, ni-nBi-1, by=1)
      grPi = grPi - sum(1/(1 - pi + k*psi))
    }
    
    dl.dpi[i] = grPi
  }

  # ---------------------------------------------
  # d_logL/d_log_kappa
  # ---------------------------------------------

  dzeta.dkappa = (gamma - eta)*(1 - rhos)*rhos/(delta^2)
  
  dl.dkappa    = kappa*sum(dl.dpi * dpi.dzeta * dzeta.dkappa)

  # ---------------------------------------------
  # d_logL/d_log_eta
  # ---------------------------------------------

  if(H0=="eta"){
    dl.deta = 0.0
  }else{
    dzeta.deta = cs
    dl.deta    = eta*sum(dl.dpi * dpi.dzeta * dzeta.deta)
  }

  # ---------------------------------------------
  # d_logL/d_log_gamma
  # ---------------------------------------------

  if(H0=="gamma"){
    dl.dgamma = 0.0
  }else{
    dzeta.dgamma = 1 - cs
    dl.dgamma    = gamma*sum(dl.dpi * dpi.dzeta * dzeta.dgamma)
  }

  if(H0==""){
    grd = c(dl.dkappa, dl.deta, dl.dgamma)
  }else if(H0=="eta"){
    grd = c(dl.dkappa, dl.dgamma)
  }else if(H0=="gamma"){
    grd = c(dl.dkappa, dl.deta)
  }

  #-----------------------------------------------------------
  # calculate Hessian matrix of ASE likelihood with respect
  # to log(kappa), log(eta), and log(gamma)
  #-----------------------------------------------------------

  if(Hessian){
    # ---------------------------------------------
    # d2l.dpi2
    # ---------------------------------------------
    
    d2l.dpi2 = rep(0, length(nB))
    
    for(i in 1:length(nB)){
      
      grPi = 0
      
      ni   = nTotal[i]
      nBi  = nB[i]
      pi   = pis[i]
      
      if(nBi > 0){
        k    = seq(0, nBi-1, by=1)
        grPi = grPi - sum(1/((pi + k*psi)^2))
      }
      
      if(nBi < ni){
        k    = seq(0, ni-nBi-1, by=1)
        grPi = grPi - sum(1/((1 - pi + k*psi)^2))
      }
      
      d2l.dpi2[i] = grPi
    }
    
    # ---------------------------------------------
    # other terms
    # ---------------------------------------------

    d2pi.dzeta2     = -2.0 / ((1 + zeta)^3)
    d2zeta.dkappa2  = -2.0 * dzeta.dkappa * rhos / delta
    dcs.dkappa      = -(1 - rhos)*rhos/(delta^2)

    dpi.dzeta.sq    = dpi.dzeta * dpi.dzeta
    dzeta.dkappa.sq = dzeta.dkappa * dzeta.dkappa

    d2lxdpi.sq = d2l.dpi2 * dpi.dzeta.sq
    dlxd2pi    = dl.dpi * d2pi.dzeta2
    dlxdpi     = dl.dpi * dpi.dzeta
    vec.d2l.dl = d2lxdpi.sq + dlxd2pi

    sum1 = sum(d2lxdpi.sq * dzeta.dkappa.sq)
    sum2 = sum(dlxd2pi * dzeta.dkappa.sq)
    sum3 = sum(dlxdpi * d2zeta.dkappa2)

    d2l.dkappa2 = dl.dkappa + kappa * kappa * (sum1 + sum2 + sum3)
    
    if(H0!="eta"){
      d2l.deta2 = dl.deta + eta * eta * sum(cs*cs*vec.d2l.dl)
      
      sum1 = sum(vec.d2l.dl * dzeta.dkappa * cs)
      sum2 = sum(dlxdpi * dcs.dkappa)
      
      d2l.dkappa.deta = kappa * eta * (sum1 + sum2)

    }
    
    if(H0!="gamma"){
      
      d2l.dgamma2 = dl.dgamma + gamma * gamma * sum((1-cs)*(1-cs)*vec.d2l.dl)
      
      sum1 = sum(vec.d2l.dl * dzeta.dkappa * (1 - cs))
      sum2 = sum(dlxdpi * dcs.dkappa)
      
      d2l.dkappa.dgamma = kappa * gamma * (sum1 - sum2)

    }
    
    if(H0==""){
      d2l.deta.dgamma = eta * gamma * sum(cs * (1-cs) * vec.d2l.dl)
      
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
# Hessian of ASE likelihood with respect to
# log(kappa), log(eta), and log(gamma)
# -----------------------------------------------------------------

neg.HessianASE.TS <- function(para, nB, nTotal, psi, rhos, H0){
  
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
  pis    = zeta/(1 + zeta)
  
  dpi.dzeta = 1/((1 + zeta)^2)
  
  # ---------------------------------------------
  # dl.dpi
  # ---------------------------------------------
  
  dl.dpi = rep(0, length(nB))
  
  for(i in 1:length(nB)){
    
    grPi = 0
    
    ni   = nTotal[i]
    nBi  = nB[i]
    pi   = pis[i]
    
    if(nBi > 0){
      k    = seq(0, nBi-1, by=1)
      grPi = grPi + sum(1/(pi + k*psi))
    }
    
    if(nBi < ni){
      k    = seq(0, ni-nBi-1, by=1)
      grPi = grPi - sum(1/(1 - pi + k*psi))
    }
    
    dl.dpi[i] = grPi
  }
  
  # ---------------------------------------------
  # d_logL/d_log_kappa
  # ---------------------------------------------
  
  dzeta.dkappa = (gamma - eta)*(1 - rhos)*rhos/(delta^2)
  
  dl.dkappa    = kappa*sum(dl.dpi * dpi.dzeta * dzeta.dkappa)
  
  # ---------------------------------------------
  # d_logL/d_log_eta
  # ---------------------------------------------
  
  if(H0=="eta"){
    dl.deta = 0.0
  }else{
    dzeta.deta = cs
    dl.deta    = eta*sum(dl.dpi * dpi.dzeta * dzeta.deta)
  }
  
  # ---------------------------------------------
  # d_logL/d_log_gamma
  # ---------------------------------------------
  
  if(H0=="gamma"){
    dl.dgamma = 0.0
  }else{
    dzeta.dgamma = 1 - cs
    dl.dgamma    = gamma*sum(dl.dpi * dpi.dzeta * dzeta.dgamma)
  }
  
  #-----------------------------------------------------------
  # calculate Hessian matrix of ASE likelihood with respect
  # to log(kappa), log(eta), and log(gamma)
  #-----------------------------------------------------------
  
  # ---------------------------------------------
  # d2l.dpi2
  # ---------------------------------------------
  
  d2l.dpi2 = rep(0, length(nB))
  
  for(i in 1:length(nB)){
    
    grPi = 0
    
    ni   = nTotal[i]
    nBi  = nB[i]
    pi   = pis[i]
    
    if(nBi > 0){
      k    = seq(0, nBi-1, by=1)
      grPi = grPi - sum(1/((pi + k*psi)^2))
    }
    
    if(nBi < ni){
      k    = seq(0, ni-nBi-1, by=1)
      grPi = grPi - sum(1/((1 - pi + k*psi)^2))
    }
    
    d2l.dpi2[i] = grPi
  }
  
  # ---------------------------------------------
  # other terms
  # ---------------------------------------------
  
  d2pi.dzeta2     = -2.0 / ((1 + zeta)^3)
  d2zeta.dkappa2  = -2.0 * dzeta.dkappa * rhos / delta
  dcs.dkappa      = -(1 - rhos)*rhos/(delta^2)
  
  dpi.dzeta.sq    = dpi.dzeta * dpi.dzeta
  dzeta.dkappa.sq = dzeta.dkappa * dzeta.dkappa
  
  d2lxdpi.sq = d2l.dpi2 * dpi.dzeta.sq
  dlxd2pi    = dl.dpi * d2pi.dzeta2
  dlxdpi     = dl.dpi * dpi.dzeta
  vec.d2l.dl = d2lxdpi.sq + dlxd2pi
  
  sum1 = sum(d2lxdpi.sq * dzeta.dkappa.sq)
  sum2 = sum(dlxd2pi * dzeta.dkappa.sq)
  sum3 = sum(dlxdpi * d2zeta.dkappa2)
  
  d2l.dkappa2 = dl.dkappa + kappa * kappa * (sum1 + sum2 + sum3)
  
  if(H0!="eta"){
    d2l.deta2 = dl.deta + eta * eta * sum(cs*cs*vec.d2l.dl)
    
    sum1 = sum(vec.d2l.dl * dzeta.dkappa * cs)
    sum2 = sum(dlxdpi * dcs.dkappa)
    
    d2l.dkappa.deta = kappa * eta * (sum1 + sum2)
    
  }
  
  if(H0!="gamma"){
    
    d2l.dgamma2 = dl.dgamma + gamma * gamma * sum((1-cs)*(1-cs)*vec.d2l.dl)
    
    sum1 = sum(vec.d2l.dl * dzeta.dkappa * (1 - cs))
    sum2 = sum(dlxdpi * dcs.dkappa)
    
    d2l.dkappa.dgamma = kappa * gamma * (sum1 - sum2)
    
  }
  
  if(H0==""){
    d2l.deta.dgamma = eta * gamma * sum(cs * (1-cs) * vec.d2l.dl)
    
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

# neg.gradASE.TS(para, nB, nTotal, psi, rhos, H0, Hessian=FALSE)


plotItASE <- function(nB, nTotal, psi, rhos, kappa, eta, gamma){
  
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
      ll[i]  = -neg.logLASE.TS(para, nB, nTotal, psi, rhos, H0="")
      gs[i,] = -neg.gradASE.TS(para, nB, nTotal, psi, rhos, H0="")
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
# ASE model
#-------------------------------------------------------------

aseR.TS <- function(y1, y2, z, rhos, H0, para0=NULL, nIter=100,
                    traceIt=FALSE, method=c("BFGS", "Newton", "nlminb"),
                    piEst=FALSE, convergence=1e-4,
                    min.ASE.Total=8, min.nASE=10)
{
  
  if(any(y1 != as.integer(y1))){
    stop("y1 must be integers\n")
  }
  
  if(any(y2 != as.integer(y2))){
    stop("y2 must be integers\n")
  }
  
  if(any(y1 < 0)){
    stop("y1 must non-negative\n")
  }
  
  if(any(y2 < 0)){
    stop("y2 must non-negative\n")
  }

  if(! all(unique(z) %in% c(0,1,3,4))){
    stop("z has invlaid values\n")
  }
  
  nn1 = length(y1)
  if( length(y2) != nn1 || length(z) != nn1 || length(rhos) != nn1){
    stop("y1, y2, z, and rhos must have the same length\n")
  }
  
  #----------------------------------------------------------
  # check the number of allele-specific reads
  #----------------------------------------------------------

  # for ASE data
  w2use = which((y1 + y2) >= min.ASE.Total)
  
  if(length(w2use) < min.nASE){
    stop("no enough samples for ASE analysis\n")
  }

  y1   = y1[w2use]
  y2   = y2[w2use]
  z    = z[w2use]
  rhos = rhos[w2use]
  
  #----------------------------------------------------------
  # problem setup
  #----------------------------------------------------------

  nTotal   = y1 + y2
  nB       = y2
  nB[z==3] = y1[z==3]
  
  # for ASE data
  wAA = which(z==0)
  wAB = which(z==1 | z==3)
  wBB = which(z==4)
  
  wHz = which(z==0 | z==4)
  pisHz = rep(0.5, length(wHz))
  
  #----------------------------------------------------------
  # initial estimation of psi
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
    
  }else{
    
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
  
  #-----------------------------------------------------------
  # initial estimate of psi
  #-----------------------------------------------------------

  cs   = (1-rhos)/(1-rhos + rhos*kappa)
  zeta = cs*eta + (1-cs)*gamma
  pis = zeta/(zeta + 1)
  pis[which(z==0 | z==4)] = 0.5

  f.lower = gradASE.psi(1e-5, nB, nTotal, pis)
  f.upper = gradASE.psi(100,  nB, nTotal, pis)
  
  u1 = uniroot(gradASE.psi, interval=c(1e-5, 100), nB, nTotal, pis,
                f.lower=f.lower, f.upper=f.upper)

  psi0 = u1$root
  
  #-----------------------------------------------------------
  # initilize likelihood as -Inf
  #-----------------------------------------------------------

  logLik  = NULL
  logLik0 = -Inf

  if(traceIt)
  message(sprintf("Initial: kappa=%.2e, eta=%.2e, gamma=%.2e, psi0=%.2e\n",
  kappa, eta, gamma, psi0))
  
  #-----------------------------------------------------------
  # take the set of samples with genotype AB
  #-----------------------------------------------------------
  
  nB.AB     = nB[wAB]
  nTotal.AB = nTotal[wAB]
  rhos.AB   = rhos[wAB]

  # neg.gradASE.TS(para, nB, nTotal, psi, rhos, H0, Hessian=FALSE)

  for(g in 1:nIter){
    
    parDiff = 0.0

    #--------------------------------------------------------
    # estimate parametrs first
    #--------------------------------------------------------
    
    if(method=="BFGS"){
      
      o1 = optim(log(para0), fn=neg.logLASE.TS, gr=neg.gradASE.TS, nB=nB.AB,
      nTotal=nTotal.AB, psi=psi0, rhos=rhos.AB, H0=H0,
      method="BFGS", control=list(maxit=10, trace=0))
      
      # o1$convergence = 0 means success
      # and o1$convergence = 1 means reaching iteration limit
      if(o1$convergence != 0 && o1$convergence != 1){
        warning("g =", g, " fail BFGS\n")
      }
      
      para1     = exp(o1$par)
      logLik1AB = -o1$value
      
    }else if(method=="Newton"){
      
      n1 = nlm(f=neg.logLASE.TS, p=log(para0), nB=nB.AB,
                nTotal=nTotal.AB, psi=psi0, rhos=rhos.AB,
                H0=H0, gH=TRUE, iterlim=10)
      # o1$convergence = 0 means success
      # and o1$convergence = 1 means reaching iteration limit
      if(n1$code == 3 || n1$code == 5){
        warning("g =", g, " nlm may not converge.\n")
      }
      
      para1     = exp(n1$estimate)
      names(para1) = names(para0)
      logLik1AB = -n1$minimum
      
    }else if(method=="nlminb"){
      
      nl1 = nlminb(log(para0), neg.logLASE.TS, gradient = neg.gradASE.TS,
                  hessian=neg.HessianASE.TS, nB=nB.AB,
                  nTotal=nTotal.AB, psi=psi0, rhos=rhos.AB,
                  H0=H0, control = list(iter.max=10))
      
      para1     = exp(nl1$par)
      logLik1AB = -nl1$objective
      
    }else{
      stop("invalid value for method\n.")
      
    }
    
    # loglikBB(nB, nTotal, pis, psi)
    logLik1Hz = loglikBB(nB[wHz], nTotal[wHz], pisHz, psi0)
    logLik1   = logLik1Hz + logLik1AB
    
    if(traceIt){
      message(sprintf("\ng=%d", g))
      message(sprintf("  Estimate para: logLik_ASE=(%.3e, %.3e)",
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
    # estimate psi
    #-----------------------------------------------------------
    
    cs   = (1-rhos)/(1-rhos + rhos*kappa)
    zeta = cs*eta + (1-cs)*gamma
    pis  = zeta/(zeta + 1)
    pis[which(z==0 | z==4)] = 0.5
    
    f.lower = gradASE.psi(1e-5, nB, nTotal, pis)
    f.upper = gradASE.psi(100,  nB, nTotal, pis)
    
    u1 = uniroot(gradASE.psi, interval=c(1e-5, 100), nB, nTotal, pis,
    f.lower=f.lower, f.upper=f.upper)
    
    psi1 = u1$root

    logLik1 = loglikBB(nB, nTotal, pis, psi1)
    
    if((logLik1 - logLik0)/abs(logLik0) < -1e-5){
      stop("liklihood decreases for estimating psi\n")
    }
    
    if(parDiff < abs(psi1 - psi0)){
      parDiff = abs(psi1 - psi1)
    }
    
    if(traceIt){
      message(sprintf("Estimate psi: psi=%.2e, logLikTReC=(%.3e, %.3e)",
      psi1, logLik0, logLik1))
    }
    
    psi0 = psi1
    logLik0 = logLik1
    
    #-----------------------------------------------------------
    # check convergence
    #-----------------------------------------------------------
    
    logLik = c(logLik, logLik1)
    
    if(traceIt){
      str = "g=%d, parDiff=%.2e, kappa=%.2f, eta=%.2f, gamma=%.2f, psi0=%.2e"
      message(sprintf(str, g, parDiff, kappa, eta, gamma, psi0))
    }
    
    if(parDiff < convergence) break
    
  }
  
  l1 = list(para=para0, psi=psi0, logLik=logLik)
  
  if(piEst){
    l1$piEst = pis
  }

  l1
}

#-------------------------------------------------------------
# ASE model
# z1 is genotype data coded as 0, 1, and 2
#-------------------------------------------------------------

# aseR.TS (y1, y2, z, rhos, H0, para0=NULL, nIter=100,
#         traceIt=FALSE, method=c("BFGS", "Newton", "nlminb"),
#         piEst=FALSE, convergence=1e-4)

aseR.TS.test <- function(y1, y2, z, rhos, nIter=100,
                        method=c("BFGS", "Newton", "nlminb"),
                        traceIt=FALSE)
{
  
  t1      = aseR.TS(y1, y2, z, rhos, nIter=nIter, H0="",
                    method=method, traceIt=traceIt)
  
  t0eta   = aseR.TS(y1, y2, z, rhos, nIter=nIter, H0="eta",
                    para0=t1$para[c(1,3)], method=method,
                    traceIt=traceIt)
  
  t0gamma = aseR.TS(y1, y2, z, rhos, nIter=nIter, H0="gamma",
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
  
  paraEst = list(para=t1$para, phi=t1$phi, bs=t1$bs)
  
  list(mle=paraEst, pvalue=pvs)
}


