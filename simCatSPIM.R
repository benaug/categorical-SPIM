e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}
simCatSPIM <-
  function(N=120,lam0=NA,p0=NA,sigma=0.50,theta=NA,K=10,X=X,buff=3,obstype="poisson",n.cat=n.cat,
           pID=pID,gamma=gamma,IDcovs=IDcovs,seed=NA){
    
    if(!is.na(seed)){
      set.seed(seed)
    }
    
    # simulate a population of activity centers
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    J<- nrow(X)
    
    #simulate IDcovs
    G.true=matrix(NA,nrow=N,ncol=n.cat) #all IDcovs in population.
    for(i in 1:N){
      for(j in 1:n.cat){
        G.true[i,j]=sample(IDcovs[[j]],1,prob=gamma[[j]])
      }
    }
    if(dim(unique(G.true))[1]!=N){
      print(paste("simulated",
                  length(unique(G.true[duplicated(G.true)]))+sum(duplicated(G.true)),"duplicated IDcovs in population"))
    }
    # Capture individuals
    y=array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      if(is.na(p0))stop("must supply p0 for bernoulli obsmod")
      pd=p0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rbinom(1,1,pd[i,j])
          }
        }
      }
    }else if(obstype=="poisson"){
      if(is.na(lam0))stop("must supply p0 for poisson obsmod")
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rpois(1,lamd[i,j])
          }
        }
      } 
    }else if(obstype=="negbin"){
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rnbinom(1,mu=lamd[i,j],size=theta)
          }
        }
      } 
    }else{
      stop("obstype not recognized")
    }

    #discard uncaptured inds and aggregate true IDcovs for all samples, keeping track of where they came from with A matrix (used to put observed data back together)
    caught=which(apply(y,c(1),sum)>0)
    y.true=y
    y=y[caught,,]
    if(K==1){
      y=array(y,dim=c(dim(y),1))
    }
    n=length(caught)
    n.samples=sum(y)
    G.cap=matrix(NA,nrow=n.samples,ncol=n.cat)
    idx=1
    A=array(0,dim=c(dim(y),n.samples))  #configuration matrix: indicator matrix for which individual i occassion j  trap k corresponds to sample l. used to convert corrupt IDcovs to corrupt capture history
    for(i in 1:length(caught)){ #loop through inds (uncaptured already removed)
      for(j in 1:J){ #then traps
        for(k in 1:K){ #then occasions
          if(y[i,j,k]>0){ #is there at least one sample here?
            for(l in 1:y[i,j,k]){ #then samples
              G.cap[idx,]=G.true[caught[i],]
              A[i,j,k,idx]=1
              idx=idx+1
            }
          }
        }
      }
    }
    ycap=aperm(apply(A,c(2,3,4),sum),c(3,1,2))
    #Amplification failure
    G.drop=G.cap #n.samples x n.cat
    for(j in 1:n.cat){
      G.drop[which(rbinom(n.samples,1,pID[j])==0),j]=0 #0 is dropout
    }
    G.obs=unique(G.drop)
    nobs=nrow(G.obs)
    yobs=array(0,dim=c(nobs,J,K))
    ID=rep(NA,n.samples)
    for(i in 1:n.samples){
      for(j in 1:nobs){
        if(all(G.drop[i,]==G.obs[j,])){
          ID[i]=j
          next
        }
      }
    }
    yobs=array(0,dim=c(max(ID),J,K))
    for(i in 1:n.samples){
      map=which(A[,,,i]==1,arr.ind=TRUE)
      yobs[ID[i],map[2],map[3]]=yobs[ID[i],map[2],map[3]]+1
    }
    
    ycheck=array(0,dim=dim(y))
    for(i in 1:n.samples){
      idx2=which(A[,,,i]>0,arr.ind=TRUE)
      ycheck[idx2[1],,]=ycheck[idx2[1],,]+A[idx2[1],,,i]
    }
    if(!all(ycheck==y)){
      stop("not all y==ycheck")
    }
    #Make ID constraint matrix
    # n.samples=sum(y.obs)
    constraints=matrix(1,nrow=n.samples,ncol=n.samples)
    for(i in 1:n.samples){
      for(j in 1:n.samples){
        guys1=which(G.drop[i,]!=0)
        guys2=which(G.drop[j,]!=0)
        comp=guys1[which(guys1%in%guys2)]
        if(any(G.drop[i,comp]!=G.drop[j,comp])){
          constraints[i,j]=0
        }
      }
    }
    #check constraints
    a=which(constraints==1,arr.ind=TRUE)#consistent
    for(i in 1:nrow(a)){
      comp=G.drop[a[i,1],]>0&G.drop[a[i,2],]>0
      if(!all(G.drop[a[i,1],comp]==G.drop[a[i,2],comp])){
        stop("Error in constraint matrix")
      }
    }
    a=which(constraints==0,arr.ind=TRUE)#inconsistent
    if(length(a)>1){
      for(i in 1:nrow(a)){
        comp=G.drop[a[i,1],]>0&G.drop[a[i,2],]>0
        if(all(G.drop[a[i,1],comp]==G.drop[a[i,2],comp])){
          stop("Error in constraint matrix")
        }
      }
    }
    A=apply(A,c(1,4),sum)
    IDtrue=rep(NA,n.samples)
    for(i in 1:n.samples){
      IDtrue[i]=which(A[,i]==1)
    }
    
    #observed capture data can be represented by site and occasion of each count member
    #different than how observed data is represented in catSPIM paper
    this.j=this.k=rep(NA,n.samples)
    for(i in 1:n.samples){
      tmp=which(ycap[i,,]==1,arr.ind=TRUE)
      this.j[i]=tmp[1]
      this.k[i]=tmp[2]
    }
    out<-list(y=y,y.obs=ycap,y.true=y.true,G.true=G.true,G.cap=G.cap,G.obs=G.drop,IDlist=list(n.cat=n.cat,IDcovs=IDcovs),
              IDtrue=IDtrue,X=X,K=K,buff=buff,constraints=constraints,obstype=obstype,s=s,n=nrow(y),
              this.j=this.j,this.k=this.k,K=K,seed=seed)
    return(out)
  }