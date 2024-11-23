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
      return(0)
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
    capcounts <- numeric(M, value = 0)
    for(i in 1:M){
      capcounts[i] <- sum(y.true[i,1:J])
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
    IDups <- control$IDups
    calcNodes <- model$getDependencies("y.true")
  },
  
  run = function() {
    z <- model$z
    G.true <- model$G.true
    for(up in 1:IDups){
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
          rem <- which(!possible)
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
GSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    M <- control$M
    n.cat <- control$n.cat
    n.levels <- control$n.levels
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    for(l in 1:n.cat){
      swap <- which(model$G.latent[1:M,l]==1)
      for(i in 1:length(swap)){
        model$G.true[swap[i],l] <<- rcat(1,model$gammaMat[l,1:n.levels[l]])
      }
    }
    # update logProb
    model_lp_proposed <- model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    z.ups <- control$z.ups
    y.nodes <- control$y.nodes
    lam.nodes <- control$lam.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected call
      if(updown==0){#subtract
        # find all z's currently on
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        #prereject turning off individuals currently allocated samples
        if(model$capcounts[pick]>0){#is this an individual with samples?
          reject <- TRUE
        }
        if(!reject){
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])
          
          #propose new N/z
          model$N[1] <<-  model$N[1] - 1
          model$z[pick] <<- 0
          
          #turn lam off
          model$calculate(lam.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn lam on
          model$calculate(lam.nodes[pick])
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)