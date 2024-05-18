e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.data.catSPIM <- function(data=NA,M=NA,inits=inits,obstype="poisson"){
  library(abind)
  this.j <- data$this.j
  this.k <- data$this.k
  X <- as.matrix(data$X)
  J <- nrow(X)
  K <- data$K
  G.obs <- data$G.obs
  n.samples <- length(this.j)
  if(is.matrix(G.obs)){
    n.cat <- ncol(G.obs)
  }else{
    n.cat <- 1
  }
  n.levels <- unlist(lapply(data$IDlist$IDcovs,length))
  #state space extent
  buff <- data$buff
  xlim <- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim <- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  #make constraints due to ID covs to initialize data
  constraints <- matrix(1,nrow=n.samples,ncol=n.samples)
  for(i in 1:n.samples){
    for(j in 1:n.samples){
      guys1 <- which(G.obs[i,]!=0)
      guys2 <- which(G.obs[j,]!=0)
      comp <- guys1[which(guys1%in%guys2)]
      if(any(G.obs[i,comp]!=G.obs[j,comp])){
        constraints[i,j] <- 0
      }
    }
  }
  #If bernoulli data, add constraints that prevent y.true[i,j,k]>1
  binconstraints <- FALSE
  if(obstype=="bernoulli"){
    for(i in 1:n.samples){
      for(j in 1:n.samples){
        if(i!=j){
          if(this.j[i]==this.j[j]&this.k[i]==this.k[j]){
            constraints[i,j] <- 0 #can't combine samples from same trap and occasion in binomial model
            constraints[j,i] <- 0
            binconstraints <- TRUE
          }
        }
      }
    }
  }
  
  #Build y.true
  y.true <- array(0,dim=c(M,J,K))
  ID <- rep(NA,n.samples)
  idx <- 1
  for(i in 1:n.samples){
    if(idx>M){
      stop("Need to raise M to initialize y.true")
    }
    y.true2D <- apply(y.true,c(1,2),sum)
    cand <- which(y.true2D[,this.j]>0)#guys caught at same traps
    if(length(cand)>0){
      if(length(cand)>1){#if more than 1 ID to match to, choose first one
        cand <- cand[1]
      }
      #Check constraint matrix
      cands <- which(ID%in%cand)#everyone assigned this ID
      if(all(constraints[i,cands]==1)){#focal consistent with all partials already assigned
        y.true[cand,this.j[i],this.k[i]] <- y.true[cand,this.j[i],this.k[i]]+1
        ID[i] <- cand
      }else{#focal not consistent
        y.true[idx,this.j[i],this.k[i]] <- 1
        ID[i] <- idx
        idx <- idx+1
      }
    }else{#no assigned samples at this trap
      y.true[idx,this.j[i],this.k[i]] <- 1
      ID[i] <- idx
      idx <- idx+1
    }
  }
  #Check assignment consistency with constraints
  checkID <- unique(ID)
  for(i in 1:length(checkID)){
    idx <- which(ID==checkID[i])
    if(!all(constraints[idx,idx]==1)){
      stop("ID initialized improperly")
    }
  }
  y.true2D <- apply(y.true,c(1,2),sum)
  known.vector <- c(rep(1,max(ID)),rep(0,M-max(ID)))
  
  #Initialize z
  z <- 1*(apply(y.true2D,1,sum)>0)
  add <- M*(0.5-sum(z)/M)
  if(add>0){
    z[sample(which(z==0),add)] <- 1 #switch some uncaptured z's to 1.
  }
  
  #Optimize starting locations given where they are trapped.
  s <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
  idx <- which(rowSums(y.true)>0) #switch for those actually caught
  for(i in idx){
    trps <- matrix(X[y.true2D[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s[i,] <- trps
    }
  }
  
  
  sigma <- inits$sigma
  D <- e2dist(s, X)
  if(obstype=="bernoulli"){
    p0 <- inits$p0
    pd <- p0*exp(-D*D/(2*sigma*sigma))
    ll.y <- dbinom(y.true2D,K,pd*z,log=TRUE)
  }else if(obstype=="poisson"){
    lam0 <- inits$lam0
    lamd <- lam0*exp(-D*D/(2*sigma*sigma))
    ll.y <- dpois(y.true2D,K*lamd*z,log=TRUE)
  }
  if(!is.finite(sum(ll.y)))stop("Starting obs model likelihood is not finite")
  
  #Initialize G.true
  G.true <- matrix(0, nrow=M,ncol=n.cat)
  for(i in 1:max(ID)){
    idx <- which(ID==i)
    if(length(idx)==1){
      G.true[i,] <- G.obs[idx,]
    }else{
      if(ncol(G.obs)>1){
        G.true[i,] <- apply(G.obs[idx,],2, max) #consensus
      }else{
        G.true[i,] <- max(G.obs[idx,])
      }
    }
  }
  G.latent <- G.true==0
  for(j in 1:n.cat){
    fix <- G.true[,j]==0
    G.true[fix,j] <- sample(1:n.levels[j],sum(fix),replace=TRUE,prob=inits$gammaMat[j,1:n.levels[j]])
  }
  return(list(y.true2D=y.true2D,y.true3D=y.true,G.true=G.true,s=s,z=z,
              ID=ID,G.latent=G.latent,n.samples=n.samples,
              xlim=xlim,ylim=ylim))
}