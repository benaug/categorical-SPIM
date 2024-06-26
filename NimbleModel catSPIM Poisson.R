NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  psi ~ dunif(0,1)
  lam0 ~ dunif(0,10)
  # sigma ~ dunif(0,100)
  sigma ~ dgamma(shape=24,rate=8)
  #category level frequencies
  for(l in 1:n.cat){
    for(k in 1:n.levels[l]){
      alpha[l,k] <- 1 #dirichlet prior parameters
    }
    gammaMat[l,1:n.levels[l]] ~ ddirch(alpha[l,1:n.levels[l]])
  }
  #--------------------------------------------------------------
  for(i in 1:M) {
    z[i] ~ dbern(psi)
    for(m in 1:n.cat){
      G.true[i,m] ~ dcat(gammaMat[m,1:n.levels[m]])
    }
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    y.true[i,1:J] ~ dPoissonVector(lam=lam[i,1:J]*K1D[1:J],z=z[i]) #vectorized obs mod
  }
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J])
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples],G.latent=G.latent[1:M,1:n.cat])
  N <- sum(z[1:M])
})#model