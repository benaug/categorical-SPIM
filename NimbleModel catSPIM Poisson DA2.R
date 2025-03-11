NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  lambda.N ~ dunif(0,1000) #expected abundance
  lam0 ~ dunif(0,10) #baseline detection rate
  sigma ~ dunif(0,10) #detection spatial scale parameter
  #category level frequencies
  for(m in 1:n.cat){
    for(k in 1:n.levels[m]){
      alpha[m,k] <- 1 #dirichlet prior parameters
    }
    gammaMat[m,1:n.levels[m]] ~ ddirch(alpha[m,1:n.levels[m]])
  }
  #--------------------------------------------------------------
  N ~ dpois(lambda.N) #realized abundance
  #data augmentation "under the hood", jointly update N/z, 
  #no distribution induced on z, just turns obsmod on/off, used in y.true/ID update
  for(i in 1:M) {
    for(m in 1:n.cat){
      G.true[i,m] ~ dcat(gammaMat[m,1:n.levels[m]])
    }
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.true[i,1:J] ~ dPoissonVector(lam=lam[i,1:J]*K1D[1:J],z=z[i]) #vectorized obs mod
  }
  #calculate number of inds captured and abundance
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M) #intermediate object
  #must use G.latent somewhere to make nimble happy. Sticking it here, not used in function.
  n <- Getncap(capcounts=capcounts[1:M],G.latent=G.latent[1:M,1:n.cat])
})#model