
#---------------------------------------------------------------
# log likelihood of TReC and ASE, for expression from
# tumor tissue with mixture of normal and tumor cells
#
# ND = no derivative, use "Nelder-Mead" method for maximization
#
# kappa = expression of tumor vs. normal for allelel A
# eta   = eQTL effect at normal cells
# gamma = eQTL effect at tumor cells
# Xs    = covariates: (intercept, read-depth, other covariates)
# bs    = coefficients for Xs
# phi   = over-dispersion parameter for TReC
# psi   = over-dispersion parameter for ASE
# rho   = tumor cell proportion
# y     = TReC
# y1    = ASE of haplotype 1
# y2    = ASE of haplotype 2
# z     = genotype 0 (AA), 1 (AB), 3 (BA) or 4(BB),
#         where A and B are two alleles of the SNP.
#         AB and BA are different since AB means the 1st/2nd
#         haplotype harbors the A/B allele, respectively,
#         while BA means the 1st/2nd haplotype harbors the
#         B/A allele, respectively.
#
# input data: rhos, y, y1, y2, Xs, z
#---------------------------------------------------------------

neg.loglik.TS.ND <- function(paras, y, nB, nTotal, Xs, z, zAE, rhos,
                             rhosAE, H0)
{
  
  if(H0 == ""){
    
    paraExp  = exp(paras[1:5])
    bs       = paras[6:length(paras)]
    
    kappa = paraExp[1]
    eta   = paraExp[2]
    gamma = paraExp[3]
    phi   = paraExp[4]
    psi   = paraExp[5]
    
  }else if(H0 == "eta"){
    
    eta     = 1.0
    paraExp = exp(paras[1:4])
    bs      = paras[5:length(paras)]
    
    kappa = paraExp[1]
    gamma = paraExp[2]
    phi   = paraExp[3]
    psi   = paraExp[4]
    
  }else if(H0 == "gamma"){

    gamma   = 1.0
    paraExp = exp(paras[1:4])
    bs      = paras[5:length(paras)]
    
    kappa = paraExp[1]
    eta   = paraExp[2]
    phi   = paraExp[3]
    psi   = paraExp[4]

  }else{
    stop("invalid value for H0!\n")
  }
  
  # for ASE data
  cs     = (1-rhosAE)/(1-rhosAE + rhosAE*kappa)
  zeta   = cs*eta + (1-cs)*gamma
  piBB   = zeta/(zeta + 1)
  
  piBB[which(zAE==0 | zAE==4)] = 0.5

  logASE = loglikBB(nB, nTotal, piBB, psi)
  
  # for TReC data
  wAA = which(z==0)
  wAB = which(z==1 | z==3)
  wBB = which(z==4)
  
  muTReC = exp(Xs %*% bs)*(1 - rhos + rhos *kappa)
  
  if(length(wAB) > 0){
    muTReC[wAB] = muTReC[wAB]*(zeta[wAB]+1)/2
  }
  
  if(length(wBB) > 0){
    muTReC[wBB] = muTReC[wBB]*zeta[wBB]
  }
  
  logTReC = loglikNB(y, muTReC, phi)
  
  -(logASE + logTReC)
}

#-------------------------------------------------------------
# joint model, ND: no derivative
#-------------------------------------------------------------

trecaseR.TS.ND <- function(dat, Xs, z, rhos, paras=NULL,
                           min.ASE.Total=8, min.nASE=10)
{
  
  #-----------------------------------------------------------
  # check parameters
  #-----------------------------------------------------------

  if(!is.data.frame(dat)){
    stop("dat is not a data.frame\n")
  }
  
  if(! all(c("y", "y1", "y2") %in% names(dat))){
    stop("dat should include 'y', 'y1', and 'y2'\n")
  }

  y  = dat$y
  y1 = dat$y1
  y2 = dat$y2
  
  if(length(y) != length(z)){
    stop("y and z have different lengths\n")
  }
  
  if(length(y) != length(rhos)){
    stop("y and rhos have different lengths\n")
  }
  
  if(length(y) != length(y1)){
    stop("y and y1 have different lengths\n")
  }
  
  if(length(y1) != length(y2)){
    stop("y1 and y2 have different lengths\n")
  }
  
  if(! all(z %in% c(0,1,3,4))){
    stop("z must have values of 0, 1, 3, or 4\n")
  }
  
  #-----------------------------------------------------------
  # check number of allele-specific reads for ASE data
  #-----------------------------------------------------------

  nTotal   = y1 + y2
  nB       = y2
  nB[z==3] = y1[z==3]
  
  # for ASE data
  w2use = which(nTotal >= min.ASE.Total)
  
  if(length(w2use) < min.nASE){
    stop("no enough samples for ASE analysis\n")
  }
  
  nB     = nB[w2use]
  nTotal = nTotal[w2use]
  zAE    = z[w2use]
  rhosAE = rhos[w2use]
  
  #-----------------------------------------------------------
  # MLE
  #-----------------------------------------------------------
  
  names(paras)[1:5] = c("kappa", "eta", "gamma", "phi", "psi")
  names(paras)[6:length(paras)] = paste("b", 1:(length(paras)-6+1), sep="")

  o1 = optim(par=paras, fn=neg.loglik.TS.ND, y, nB, nTotal, Xs=Xs, z=z,
             zAE=zAE, rhos=rhos, rhosAE=rhosAE, H0="")

  paras.eta   = paras[-2]
  paras.gamma = paras[-3]

  o1.eta   = optim(par=paras.eta, fn=neg.loglik.TS.ND, y, nB, nTotal, Xs=Xs,
                   z=z, zAE=zAE, rhos=rhos, rhosAE=rhosAE, H0="eta")

  o1.gamma = optim(par=paras.gamma, fn=neg.loglik.TS.ND, y, nB, nTotal, Xs=Xs,
                   z=z, zAE=zAE, rhos=rhos, rhosAE=rhosAE, H0="gamma")

  paraEst = o1$par
  paraEst[1:5] = exp(paraEst[1:5])

  l1       = -o1$value
  l0.eta   = -o1.eta$value
  l0.gamma = -o1.gamma$value
  
  lrt.eta   = 2*(l1 - l0.eta)
  lrt.gamma = 2*(l1 - l0.gamma)

  pv.eta    = pchisq(lrt.eta, df=1, lower.tail=FALSE)
  pv.gamma  = pchisq(lrt.gamma, df=1, lower.tail=FALSE)

  pvs = data.frame(LRS=c(lrt.eta, lrt.gamma), pval=c(pv.eta, pv.gamma))
  rownames(pvs) = c("eta", "gamma")

  list(mle=paraEst, pvalue=pvs, loglik=l1)
  
}

# -----------------------------------------------------------------
# log likelihood of TReCASE data with respect to log(kappa),
# log(eta), and log(gamma)
# -----------------------------------------------------------------

neg.logLik.TS <- function(paras, y, nB, nTotal, z1, nu0, phi, psi,
                          rhos, rhosAE, H0)
{
  # for ASE data
  nLogLikASE  = neg.logLASE.TS(paras, nB, nTotal, psi, rhosAE, H0)
  
  # for TReC data
  nLogLikTReC = neg.logLTReC.TS(paras, y, z1, nu0, phi, rhos, H0)
  
  nLogLikASE + nLogLikTReC
}

# -----------------------------------------------------------------
# gradient of TReCASE log likelihood with respect to log(kappa),
# log(eta), and log(gamma)
# -----------------------------------------------------------------

neg.grad.TS <- function(paras, y, nB, nTotal, z1, nu0, phi, psi,
                        rhos, rhosAE, H0)
{
  # for ASE data
  nGradASE  = neg.gradASE.TS(paras, nB, nTotal, psi, rhosAE, H0)
  
  # for TReC data
  nGradTReC = neg.gradTReC.TS(paras, y, z1, nu0, phi, rhos, H0)

  nGradASE + nGradTReC

}

# -----------------------------------------------------------------
# Hessian of TReCASE log likelihood with respect to
# log(kappa), log(eta), and log(gamma)
# -----------------------------------------------------------------

neg.Hessian.TS <- function(paras, y, nB, nTotal, z1, nu0, phi, psi,
                           rhos, rhosAE, H0)
{
  # for ASE data
  nHessianASE  = neg.HessianASE.TS(paras, nB, nTotal, psi, rhosAE, H0)
  
  # for TReC data
  nHessianTReC = neg.HessianTReC.TS(paras, y, z1, nu0, phi, rhos, H0)
  
  nHessianASE + nHessianTReC
  
}

#-------------------------------------------------------------
# TReCASE model
#-------------------------------------------------------------

trecaseR.TS <- function(dat, Xs, z, rhos, H0, para0=NULL,
                        min.ASE.Total=8, min.nASE=10,
                        nIter=200, convergence=1e-4){
  
  #-----------------------------------------------------------
  # check parameters
  #-----------------------------------------------------------

  if(!is.data.frame(dat)){
    stop("dat is not a data.frame\n")
  }
  
  if(! all(c("y", "y1", "y2") %in% names(dat))){
    stop("dat should include 'y', 'y1', and 'y2'\n")
  }
  
  y  = dat$y
  y1 = dat$y1
  y2 = dat$y2
  
  if(length(y) != length(z)){
    stop("y and z have different lengths\n")
  }
  
  if(length(y) != length(rhos)){
    stop("y and rhos have different lengths\n")
  }
  
  if(length(y) != length(y1)){
    stop("y and y1 have different lengths\n")
  }
  
  if(length(y1) != length(y2)){
    stop("y1 and y2 have different lengths\n")
  }
  
  if(! all(z %in% c(0,1,3,4))){
    stop("z must have values of 0, 1, 3, or 4\n")
  }
  
  #-----------------------------------------------------------
  # check number of allele-specific reads for ASE data
  #-----------------------------------------------------------

  nTotal   = y1 + y2
  nB       = y2
  nB[z==3] = y1[z==3]
  
  w2use = which(nTotal >= min.ASE.Total)
  
  if(length(w2use) < min.nASE){
    stop("no enough samples for ASE analysis\n")
  }
  
  nB     = nB[w2use]
  nTotal = nTotal[w2use]
  zAE    = z[w2use]
  rhosAE = rhos[w2use]
  
  #----------------------------------------------------------
  # get initial values for the parameters
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
  # initial estimate b and phi
  #-----------------------------------------------------------
  
  z1 = z/2
  z1[z1>0.1 & z1<1.9] = 1
  
  ww1  = which(z1==1)
  ww2  = which(z1==2)
  
  gc1  = glm.control(epsilon=1e-5, maxit=25, trace = FALSE)

  cs   = (1-rhos)/(1-rhos + rhos*kappa)
  zeta = cs * eta + (1 - cs)*gamma
  
  off1 = log(1-rhos + rhos*kappa)
  off1[ww1] = off1[ww1] + log((zeta[ww1] + 1)/2)
  off1[ww2] = off1[ww2] + log(zeta[ww2])
  
  g0    = glm.nb(y ~ -1 + Xs + offset(off1),  control=gc1)
  nu0   = exp(Xs %*% g0$coef)
  phi0  = 1/g0$theta
  
  #-----------------------------------------------------------
  # initial estimate of psi
  #-----------------------------------------------------------

  zetaAE = zeta[w2use]
  pis    = zetaAE/(zetaAE + 1)
  pis[which(zAE==0 | zAE==4)] = 0.5

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

  #-----------------------------------------------------------
  # take the set of samples with genotype AB
  #-----------------------------------------------------------

  wAA = which(zAE==0)
  wAB = which(zAE==1 | zAE==3)
  wBB = which(zAE==4)

  wHz = which(zAE==0 | zAE==4)
  pisHz = rep(0.5, length(wHz))

  nB.AB     = nB[wAB]
  nTotal.AB = nTotal[wAB]
  rhos.AB   = rhos[wAB]

  #-----------------------------------------------------------
  # start iterations
  #-----------------------------------------------------------

  for(g in 1:nIter){
    
    parDiff = 0.0

    #--------------------------------------------------------
    # estimate parametrs first
    #--------------------------------------------------------
    
    # neg.logLik.TS(paras, y, nB, nTotal, z1, nu0, phi, psi,
    # rhos, rhosAE, H0)

    nl1 = nlminb(log(para0), neg.logLik.TS, gradient = neg.grad.TS,
                 hessian=neg.Hessian.TS, y=y, nB=nB.AB,
                 nTotal=nTotal.AB, z1=z1, nu0=nu0, phi=phi0,
                 psi=psi0, rhos=rhos, rhosAE=rhos.AB, H0=H0,
                 control = list(iter.max=10))
    
    para1   = exp(nl1$par)
    logLik1 = -nl1$objective
    
    # loglikBB(nB, nTotal, pis, psi)
    logLik1Hz = loglikBB(nB[wHz], nTotal[wHz], pisHz, psi0)
    logLik1   = logLik1 + logLik1Hz
    
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
    logLNB1 = loglikNB(y, nu1off, phi1)
    
    #-----------------------------------------------------------
    # estimate psi
    #-----------------------------------------------------------
    
    zetaAE = zeta[w2use]
    pis    = zetaAE/(zetaAE + 1)
    pis[which(zAE==0 | zAE==4)] = 0.5
    
    f.lower = gradASE.psi(1e-5, nB, nTotal, pis)
    f.upper = gradASE.psi(100,  nB, nTotal, pis)
    
    u1 = uniroot(gradASE.psi, interval=c(1e-5, 100), nB, nTotal, pis,
    f.lower=f.lower, f.upper=f.upper)
    
    psi1 = u1$root
    
    logLBB1 = loglikBB(nB, nTotal, pis, psi1)

    #-----------------------------------------------------------
    # check log likelihood
    #-----------------------------------------------------------

    logLik1 = logLNB1 + logLBB1
    
    if((logLik1 - logLik0)/abs(logLik0) < -1e-5){
      stop("liklihood decreases for estimating b and phi\n")
    }
    
    if(parDiff < abs(phi1 - phi0)){
      parDiff = abs(phi1 - phi0)
    }
    
    phi0 = phi1
    nu0  = exp(Xs %*% g1$coef)
    psi0 = psi1

    logLik0 = logLik1
    
    #-----------------------------------------------------------
    # check convergence
    #-----------------------------------------------------------
    
    logLik = c(logLik, logLik1)
    
    if(parDiff < convergence) break
    
  }
  
  l1 = list(para=para0, phi=phi0, psi=psi0, bs=g1$coef, logLik=logLik)
  
  l1
}

#-------------------------------------------------------------
# TReCASE model
# z1 is genotype data coded as 0, 1, and 2
#-------------------------------------------------------------


trecaseR.TS.test <- function(dat, Xs, z, rhos, nIter=200)
{
  
  
  t0eta   = trecaseR.TS(dat, Xs, z, rhos, nIter=nIter, H0="eta")
  
  t0gamma = trecaseR.TS(dat, Xs, z, rhos, nIter=nIter, H0="gamma")
  
  t1eta   = trecaseR.TS(dat, Xs, z, rhos, nIter=nIter, H0="",
                        para0=c(t0eta$para[1], 1, t0eta$para[2]))
                        
  t1gamma = trecaseR.TS(dat, Xs, z, rhos, nIter=nIter, H0="",
                        para0=c(t0eta$para, 1))

  l0eta   = t0eta$logLik[length(t0eta$logLik)]
  l0gamma = t0gamma$logLik[length(t0gamma$logLik)]
  
  l1eta   = t1eta$logLik[length(t1eta$logLik)]
  l1gamma = t1gamma$logLik[length(t1gamma$logLik)]
  
  if(l1eta > l1gamma){
    t1 = t1eta
  }else{
    t1 = t1gamma
  }
  
  l1      = t1$logLik[length(t1$logLik)]
  
  lrt.eta   = 2*(l1 - l0eta)
  lrt.gamma = 2*(l1 - l0gamma)
  
  pv.eta    = pchisq(lrt.eta, df=1, lower.tail=FALSE)
  pv.gamma  = pchisq(lrt.gamma, df=1, lower.tail=FALSE)
  
  pvs = data.frame(LRS=c(lrt.eta, lrt.gamma), pval=c(pv.eta, pv.gamma))
  rownames(pvs) = c("eta", "gamma")
  
  list(para=t1$para, phi=t1$phi, bs=t1$bs, pvalue=pvs)
}

