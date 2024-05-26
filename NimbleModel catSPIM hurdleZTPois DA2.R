NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  lambda.N ~ dunif(0,1000) #expected abundance
  logit(p0) ~ dlogis(0,1) #baseline detection probability on logit scale
  lambda ~ dunif(0,10) #parameter for counts|detection
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
    pd[i,1:J] <- GetDetectionProb(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, p0=p0, z=z[i])
    y.true[i,1:J,1:K] ~ dHurdleZTPoisMatrix(pd=pd[i,1:J],K2D=K2D[1:J,1:K],z=z[i],lambda=lambda) #vectorized obs mod
  }
  capcounts[1:M] <- Getcapcounts(y.true=y.true[1:M,1:J,1:K])
  n <- Getncap(capcounts=capcounts[1:M],ID=ID[1:n.samples],G.latent=G.latent[1:M,1:n.cat])
})#model