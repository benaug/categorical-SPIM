#This script uses a negative binomial observation model with mean parameter lambda and dispersion parameter theta
#The custom ID update must now use a Metropolis-Hastings update since the full conditional is not known.
#If you modify the model code to put covariates on the NB dispersion parameter, need to fix that in the custom ID update.

#Generally, you're going to need more informative partial IDs to get decent density estimates compared to Poisson obs mod
#Also, theta needs to be quite small (larger overdispersion) to detect overdispersion. Overdispersion is harder to
#detect compared to a typical negative binomial regression because individual by trap expected counts a funcition
#of activity center location, which also must be estimated

library(nimble)
library(coda)
source("simCatSPIM.R")
source("init.data.CatSPIM.R")
source("NimbleModel catSPIM NegBin.R")
source("NimbleFunctions catSPIM NegBin.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#simulate some data
N <- 50
lam0 <- 0.5
theta <- 0.1 #lower is more overdispersion
sigma <- 0.50
K <- 10
buff <- 3 #state space buffer. Should be at least 3 sigma.
X <- expand.grid(3:11,3:11)
obstype <- "negbin"

#categorical identity covariate stuff
n.cat <- 3  #number of ID covariates
gamma <- vector("list",n.cat) #population frequencies of each category level

#Or this for generic ID covariates
gamma <- vector("list",n.cat)
n.levels <- c(7,7,7) #Number of levels per ID covariate
for(i in 1:n.cat){
  gamma[[i]] <- rep(1/n.levels[i],n.levels[i]) #generating all equal category level frequencies
}
IDcovs <- vector("list",n.cat)
for(i in 1:length(IDcovs)){
  IDcovs[[i]] <- 1:n.levels[i]
}
pID <- rep(1,n.cat)#sample by covariate level observation probability.  e.g. loci amplification probability

#n.cat=1 with n.levels=1 will produce unmarked SCR data with no ID covariates. 
#Well, everyone has the same covariate value so they are effectively unmarked

#Simulate some data
data <- simCatSPIM(N=N,lam0=lam0,theta=theta,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype,
                n.cat=n.cat,pID=pID,gamma=gamma,
                IDcovs=IDcovs)
#look for infrequent, larger individual by trap by occasion counts indicative of overdispersion
table(data$y.true)

#What is the observed data?
#1) We have occasions and sites for each count member.
head(data$this.j)
head(data$this.k)#not used in this sampler for 2D data, but required if using 3D data
#2) We have partial ID covariates for each count member(0 is missing value)
head(data$G.obs)

#Data augmentation level
M <- 175

#trap operation matrix
J <- nrow(X)
K1D <- rep(K,J)

#set initial values for category levels
#Stuff category level probabilities into a ragged matrix.
gammaMat <- matrix(0,nrow=n.cat,ncol=max(n.levels))
for(l in 1:n.cat){
  gammaMat[l,1:n.levels[l]] <- gamma[[l]]
}
inits <- list(lam0=1,sigma=1,gammaMat=gammaMat)#initial values for lam0, sigma, and gammas to build data
nimbuild <- init.data.catSPIM(data=data,M=M,inits=inits) #watch for error building data if M needs to be larger to initialize.
G.obs.seen <- (data$G.obs!=0) #used in custom update to indicate which are observed

#inits for nimble
Niminits <- list(z=nimbuild$z,s=nimbuild$s,G.true=nimbuild$G.true,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true2D),
                 y.true=nimbuild$y.true2D,G.latent=nimbuild$G.latent,lam0=inits$lam0,sigma=inits$sigma,
                 gammaMat=gammaMat,
                 theta=1)

#constants for Nimble
constants <- list(M=M,J=J,K1D=K1D,n.samples=nimbuild$n.samples,n.cat=n.cat,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim,n.levels=n.levels)

#supply data to nimble
Nimdata <- list(y.true=matrix(NA,nrow=M,ncol=J),
              G.true=matrix(NA,nrow=M,ncol=n.cat),ID=rep(NA,nimbuild$n.samples),
              z=rep(NA,M),X=as.matrix(data$X),capcounts=rep(NA,M))

# set parameters to monitor
parameters <- c('psi','lam0','theta','sigma','N','n','gammaMat')

#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c("ID")
nt <- 1 #thinning rate
nt2 <- 50#thin more

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
#can use "nodes" argument in configureMCMC below to omit y.true and G.true that are replaced below for faster
#configuration
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt, monitors2=parameters2,thin2=nt2,useConjugacy = TRUE) 

#conf$printSamplers() #shows the samplers used for each parameter and latent variable

###Two *required* sampler replacements

##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
conf$removeSampler("y.true")
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                type = 'IDSampler',control = list(M=M,J=J,this.j=data$this.j,G.obs=data$G.obs,
                                                  G.obs.seen=G.obs.seen,n.cat=n.cat,n.samples=nrow(data$G.obs),
                                                  K1D=K1D),
                silent = TRUE)

#replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
#2 options. First one is faster, but does not allow other parameters, say lam0 and sigma, to vary
#as a function of ID covs in G.true
conf$removeSampler("G.true")
conf$addSampler(target = paste("G.true[1:",M,",1:",n.cat,"]", sep=""),
                type = 'GSampler1',
                control = list(M=M,n.cat=n.cat,n.levels=n.levels), silent = TRUE)

# #Second option allows other parameters to vary as function of G.true
# conf$removeSampler("G.true")
# for(i in 1:M){
#   for(m in 1:n.cat){ #don't need to update first cat bc it is mark status
#     conf$addSampler(target = paste("G.true[",i,",",m,"]", sep=""),
#                     type = 'GSampler2',
#                     control = list(i = i,m=m,M=M,n.cat=n.cat,n.samples=nimbuild$n.samples,
#                                    n.levels=n.levels,G.noID=nimbuild$G.noID), silent = TRUE)
#   }
# }



# ###Two *optional* sampler replacements:
# 
#replace default activity center sampler that updates x and y locations separately with a joint update
#should be a little more efficient. Slice seems better than block random walk
#Can also use "sSampler" from other repos that only tunes when z=1
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
  #                 type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
  # block RW option
  # do not adapt covariance bc samples not deterministically linked to individuals
  # longer adapt interval to average over more data configurations for each s_i
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'RW_block',control=list(adaptive=TRUE,adaptScaleOnly=TRUE,adaptInterval=500),silent = TRUE)
}

#use block update for lam0 and sigma. bc correlated posteriors.
conf$removeSampler(c("lam0","sigma"))
conf$addSampler(target = c("lam0","sigma"),
                type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)

#Might need longer adapt interval for theta if slow to converge so that it is well tuned after convergence.
#In general, it looks like tuning the theta update is difficult with this algorithm for initializing the
#latent capture data and sample IDs, which can be far from the truth.
#Providing a small initial value for theta will help.
conf$removeSampler("theta")
conf$addSampler(target = c("theta"),
                type = 'RW',control = list(adaptive=TRUE,adaptInterval=500),silent = TRUE)


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)



# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

#n is number of individuals captured. True value:
data$n
