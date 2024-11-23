#------------------------------------------------------------------
# Function for calculation detection rate
#------------------------------------------------------------------
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lam = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dpois(x, lambda = lam, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0),lam = double(1),z = double(0)) {
    returnType(double(1))
    J <- nimDim(lam)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

Getcapcounts <- nimbleFunction(
  run = function(y.true=double(2)){
    returnType(double(1))
    M <- nimDim(y.true)[1]
    J <- nimDim(y.true)[2]
    capcounts=numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i]=sum(y.true[i,1:J])
    }
    return(capcounts)
  }
)
Getncap <- nimbleFunction(
  run = function(capcounts=double(1),ID=double(1),G.latent=double(2)){ #don't need ID, but nimble requires is it used in a function 
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)

## sampler to update y[1:M,1:J] subject to G.obs constraints
IDSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    this.j <- control$this.j
    G.obs <- control$G.obs
    G.obs.seen <- control$G.obs.seen
    trapup <- control$trapup
    M <- control$M
    J <- control$J
    n.cat <- control$n.cat
    calcNodes <- model$getDependencies(target)
  },
  
  run = function() {
    z <- model$z
    G.true <- model$G.true
    for(j in 1:length(trapup)){ #only iterate through traps with samples
      lam.curr <- model$lam[,trapup[j]] #individual by trap expected counts
      lam.curr[z==0] <- 0 #can zero out z==0 here, will already be zeroed out using supplied nimble functions
      idx <- which(this.j==trapup[j]) #which samples were recorded at this trap?
      for(l in 1:length(idx)){ #update the individual identities of these samples one by one
        lam.use <- lam.curr
        #zero out lam.use for individuals whose G.true conflict with G.obs
        #this prevents assigning samples to individuals that do not match cov values
        idx2 <- which(G.obs.seen[idx[l],])#which indices of this G.obs are not missing values?
        if(length(idx2)>1){#multiple loci observed
          match <- nimLogical(M,TRUE) #start with all true, remove conflicts loci by loci
          for(l2 in 1:length(idx2)){#loop through observed loci for sample l
            match <- match[1:M]&(G.true[1:M,idx2[l2]]==G.obs[idx[l],idx2[l2]])
          }
          possible <- match
        }else if(length(idx2)==1){#single loci observed
          possible <-G.true[1:M,idx2[1]]==G.obs[idx[l],idx2[1]]
        }else{#fully latent G.obs
          possible <- nimLogical(M,TRUE) #Can match anyone with z==1
        }
        #remove individuals whose G.true conflicts with this G.obs
        rem=which(!possible)
        if(length(rem)>0){
          lam.use[rem] <- 0
        }
        #full conditional for identity update
        fullcond <- lam.use/sum(lam.use)
        #update derived parameter ID
        ID.curr <- model$ID[idx[l]]
        ID.prop <- rcat(1,fullcond)
        model$ID[idx[l]] <<- ID.prop
        #update y
        model$y.true[ID.curr,trapup[j]] <<- model$y.true[ID.curr,trapup[j]]-1 #subtract out old ID obs
        model$y.true[ID.prop,trapup[j]] <<- model$y.true[ID.prop,trapup[j]]+1 #add in new ID obs
      }
    }
    #Now update G.latent. G.latent determines which G.true can be updated later.
    #Only individuals with no samples assigned to them can have their G.true updated.
    G.true.tmp <- matrix(0, nrow=M,ncol=n.cat)
    #I only want to loop over unique(ID), but didn't spend time trying to figure out
    # and efficient way to code unique(), which is not available in NIMBLE.
    for(i in 1:length(model$ID)){
      for(l in 1:n.cat){
        if(G.obs[i,l]!=0){
          G.true.tmp[model$ID[i],l] <- G.obs[i,l]
        }
      }
    }
    model$G.latent[1:M,1:n.cat] <<-  (G.true.tmp==0)*1
    #update lp
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

## sampler to update G.true subject to constraints stemming from G.obs
#assigned to them.
#2 options here. 
#GSampler1 should be faster, but does not allow G.true to be modeled as a function of other
#parameters, e.g., sex-specific detection functions
#Just drawing from prior when a G is not fixed by an assigned sample and skipping otherwise.
GSampler1 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    M <- control$M
    n.cat <- control$n.cat
    n.levels <- control$n.levels
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    for(m in 1:n.cat){
      swap <- which(model$G.latent[1:M,m]==1)
      for(i in 1:length(swap)){
        model$G.true[swap[i],m] <<- rcat(1,model$gammaMat[m,1:n.levels[m]])
      }
    }
    # update logProb
    model_lp_proposed <- model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#Using Metropolis-Hastings here for unfixed G values so that parameters, say,
#lam0 and sigma, can vary by G.true values
GSampler2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    M <- control$M
    n.levels <- control$n.levels
    calcNodes <- model$getDependencies(target)
    i <- control$i
    m <- control$m
  },
  run = function() {
    if(model$G.latent[i,m]==1){ #skip if not latent
      model.lp.initial <- model$getLogProb(calcNodes) #initial logProb
      prop.back <- model$gammaMat[m,model$G.true[i,m]] #backwards proposal prob
      G.prop <- rcat(1,model$gammaMat[m,1:n.levels[m]])
      prop.for <- model$gammaMat[m,G.prop] #forwards proposal prob
      model$G.true[i,m] <<- G.prop #store in model
      model.lp.proposed <- model$calculate(calcNodes)#proposed logProb
      log_MH_ratio <- (model.lp.proposed+log(prop.back)) - (model.lp.initial+log(prop.for))
      
      # log_MH_ratio
      accept <- decide(log_MH_ratio)
      if(accept) {
        copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
      } else {
        copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
      }
    }
  },
  methods = list( reset = function () {} )
)