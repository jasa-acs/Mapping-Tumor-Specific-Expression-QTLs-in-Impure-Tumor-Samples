
# ----------------------------------------------------------------
# normal quantile transformation
# ----------------------------------------------------------------

normscore = function(vec) {
  len  = length(na.omit(vec))+1
  rank = rank(na.omit(vec))
  ties = (rank - floor(rank)) > 0
  new.vec = vec[!is.na(vec)] 
  new.vec[!ties]=qnorm(rank[!ties]/len)
  new.vec[ties] =0.5*(qnorm((rank[ties]+0.5)/len)+qnorm((rank[ties]-0.5)/len))
  vec[!is.na(vec)] = new.vec
  vec
}

# ----------------------------------------------------------------
# power calculation
# ----------------------------------------------------------------

powerNB <- function(n, mu, foldT, foldS, phi, psi, n.simu, alpha.h, 
                    maf, pis, rA)
{
  
  n0 = round(n*(1-maf)*(1-maf))
  n1 = round(n*maf*(1-maf))
  n3 = round(2*n*maf*(1-maf)) - n1
  n4 = n - n0 - n1 - n3
  
  z  = rep(c(0,1,3,4), times=c(n0, n1, n3, n4))
  z1 = z/2
  z1[which(z1>.1 & z1<1.9)] = 1
  
  Xs  = matrix(1, nrow=n, ncol=1)
  
  bs = log(100)
  
  eta   = 2*foldS - 1
  gamma = 2*foldT - 1
  cs   = (1-pis)/(1-pis + pis*rA)
  zeta = cs*eta + (1-cs)*gamma
  
  mu0  = exp(Xs %*% bs)
  
  muy  = mu0*(1-pis + pis*rA)
  
  wAA  = which(z==0)
  wAB  = which(z==1 | z==3)
  wBB  = which(z==4)
  
  muy[wAB] = muy[wAB] * (1 + zeta[wAB])/2
  muy[wBB] = muy[wBB] * zeta[wBB]
  
  ww0   = which(z1==0)
  ww1   = which(z1==1)
  ww2   = which(z1==2)
  
  pval1 = pval2 = pval3 = pval4 = pval5 = pval6 = pval7 = rep(1, n.simu)

  fail_ct = 0
  
  for(k in 1:n.simu){
    #-------------------------------------------------
    # Steal the seed
    #-------------------------------------------------
    curr.seed = .Random.seed
    
    message("Simulation ",k," Started")
    
    # ------------------------------------------------
    # simulate total read count
    # ------------------------------------------------
    
    y = rnegbin(n, muy, 1/phi)
    
    # ------------------------------------------------
    # simulate data: allele specific read
    # ------------------------------------------------
    
    pi    = zeta/(zeta + 1)
    pi[which(z==0 | z==4)] = 0.5
    
    alpha = pi/psi
    beta  = 1/psi - alpha
    
    nTotal = round(0.05*y)
    
    y1 = y2 = rep(NA, n)
    
    for(i in 1:n){
      if(nTotal[i] == 0){
        y1[i] = 0
        y2[i] = 0
      }else{
        y2[i] = rbetabinom.ab(1, nTotal[i], alpha[i], beta[i])
        y1[i] = nTotal[i] - y2[i]
      }
    }
    
    ww3 = which(z == 3)
    if(length(ww3) > 0){
      ytmp    = y2[ww3]
      y2[ww3] = y1[ww3]
      y1[ww3] = ytmp
    }
    
    # ------------------------------------------------
    # fit different models
    # ------------------------------------------------
    
    ## linear model
    yn = normscore(y)
    l1 = summary(lm(yn ~ z1))
    pval1[k] = l1$coef[2,4]
    
    ## linear model ala the cell type specific eQTL paper:
    l2 = summary(lm(yn ~ z1+pis+pis:z1))
    pval6[k] = l2$coef[4,4]
    
    ## Linear model ala the westra paper with different test:
    data.test = data.frame(yn=yn,z1=z1,pis=pis)
    l2.newtest= lm(yn~z1*pis,data=data.test)
    pval7[k] = summary(glht(model = l2.newtest,linfct = "z1+z1:pis=0"))$test$pvalues[1]
    
    ## trec
    t1 = trecR(y, Xs, z1, fam="negbin")
    pval2[k] = 1 - pchisq(t1$lrt,1)
    
    nTotal = y1 + y2
    wkp    = which(nTotal >= 5)
    if(length(wkp) >= 5){
      ## trecase
      z2 = z1
      z2[which(z1==2)] = 4
      t2 = trecaseR(y, y1, y2, Xs, z)
      pval3[k] = 1 - pchisq(t2$lrt,1)
    }else{
      pval3[k]  = pval3[k]
    }
    
    t4 = trecR.TS.test(y, Xs, z1, pis,method="nlminb")
    pval4[k] = t4$pvalue[2,2]
    
    currFile = sprintf("DJ_Test_FT_%s_FS_%s_psi_%s.txt", gamma, eta, psi)
    
    t5 = tryCatch(pTReCASE_multapp(Y = as.matrix(y), Y1 = as.matrix(y1),
                                   Y2 = as.matrix(y2), Z = as.matrix(z), 
                                   X = Xs, rhos = pis, F_out = currFile, 
                                   geno_pos = c(1), gene_start = c(1), 
                                   gene_end = c(1), Perform_CT_co=0),
                  error = function(e){
                    print(paste("MY_ERROR: ",e))
                    f.name = sprintf("err_foldT_%s_foldS_%s.RData",foldT,foldS)
                    save(curr.seed,file = f.name)
                    return(-10)
                    })
    if(t5==0){
      t5 = read.table(file = currFile,header = TRUE,sep = "\t", 
                      stringsAsFactors = FALSE, fill = TRUE)
      psi_idx  = which(t5$Psi>=0)
      pval5[k] = t5[psi_idx,"P_Gamma"]
      if(any(c(t5$Fail_Full,t5$Fail_Eta,t5$Fail_Gamma)%in%c(1,4))){
        fail_ct = fail_ct+1
        pval5[k] = 1
      }
    } else {
      stop("See Error Above")
    }
    
    message("Simulation ",k," Complete!")
  }
  
  pp1 = length(which(pval1 < alpha.h))/n.simu
  pp2 = length(which(pval2 < alpha.h))/n.simu
  pp3 = length(which(pval3 < alpha.h))/n.simu
  pp4 = length(which(pval4 < alpha.h))/n.simu
  pp5 = length(which(pval5 < alpha.h))/n.simu
  pp6 = length(which(pval6 < alpha.h))/n.simu
  pp7 = length(which(pval7 < alpha.h))/n.simu
  
  message("Fail Ct: ",fail_ct)
  
  c(pp1, pp2, pp3, pp4, pp5, pp6, pp7)
}

