NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  psi ~ dunif(0,1)
  lam0 ~ dunif(0,10)
  theta ~ dunif(0,25) #careful with this prior. Too much prior mass near 0 gives very strong prior weight to high overdispersion
  sigma ~ dunif(0,100)
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
    p[i,1:J] <- theta/(theta+lam[i,1:J])
    y.true[i,1:J] ~ dNBVector(p=p[i,1:J],theta=theta*K1D[1:J],z=z[i]) #vectorized obs mod. trap op: sum of NB RVs is NB with theta=N*theta
  }
  #calculate number of inds captured and abundance
  capcounts[1:M] <- Getcapcounts(ID=ID[1:n.samples],M=M) #intermediate object
  #must use G.latent somewhere to make nimble happy. Sticking it here, not used in function.
  n <- Getncap(capcounts=capcounts[1:M],G.latent=G.latent[1:M,1:n.cat])
  N <- sum(z[1:M])
})#model