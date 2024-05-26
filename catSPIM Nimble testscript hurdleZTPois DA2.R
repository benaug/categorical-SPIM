#This script uses a hurdle zero-truncated Poisson observation model
#This is a model where 1) detection is a function of distance from activity center
#and 2) given detection, counts follow a zero-truncated Poisson with a fixed lambda
#parameter that is not a function of distance from the activity center.
#We assume y.det[i,j,k] ~ Bernoulli(pd[i,j]) and
# y.count[i,j,k] ~ ZTPois(lambda*y.det[i,j,k]) (lambda is zeroed out if no detection)
#To fit the model, we marginalize over y.det.

#this version uses an alternative data augmentation approach that runs faster and allows a poisson
#prior on N.

library(nimble)
library(coda)
source("simCatSPIM.R")
source("init.data.CatSPIM.R")
source("NimbleModel catSPIM hurdleZTPois DA2.R")
source("NimbleFunctions catSPIM hurdleZTPois DA2.R")
source("sSampler.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)

#simulate some data
N <- 75
p0 <- 0.25 #baseline detection probability
lambda <- 0.25 #count parameter given detection
sigma <- 0.5
K <- 10 #number of occasions
buff <- 3 #state space buffer. Should be at least 3 sigma.
X <- expand.grid(3:11,3:11)
xlim <- range(X[,1]) + c(-buff,buff)
ylim <- range(X[,2]) + c(-buff,buff)
diff(xlim)*diff(ylim) #state space area

#categorical identity covariate stuff
n.cat <- 2  #number of ID covariates
gamma <- vector("list",n.cat)
n.levels <- c(2,5) #Number of levels per ID covariate (sex and coat)
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
data <- simCatSPIM(N=N,p0=p0,lambda=lambda,sigma=sigma,K=K,X=X,buff=buff,obstype="hurdleZTPois", #this nimble model set up for hurdleZTPois
                n.cat=n.cat,pID=pID,gamma=gamma,
                IDcovs=IDcovs)

#What is the observed data?
#1) We have occasions and sites for each count member.
head(data$this.j)
head(data$this.k)#not used in this sampler for 2D data, but required if using 3D data
#2) We have partial ID covariates for each count member(0 is missing value)
head(data$G.obs)

#Data augmentation level
M <- 250

#trap operation matrix
J <- nrow(X)
K2D <- matrix(1,J,K) #2 dimensional. Must be either 0 or 1.

#set initial values for category levels
#Stuff category level probabilities into a ragged matrix.
gammaMat <- matrix(0,nrow=n.cat,ncol=max(n.levels))
for(l in 1:n.cat){
  gammaMat[l,1:n.levels[l]] <- gamma[[l]]
}
inits <- list(p0=0.5,lambda=1,sigma=1,gammaMat=gammaMat)#ballpark initial values for lam0, sigma, and gammas to build data
nimbuild <- init.data.catSPIM(data=data,M=M,inits=inits,obstype="hurdleZTPois")
G.obs.seen <- (data$G.obs!=0) #used in custom update to indicate which are observed

#inits for nimble
Niminits <- list(z=nimbuild$z,N=sum(nimbuild$z),#must initialize N to be sum(z.init) for this data augmentation approach
                 lambda.N=sum(nimbuild$z), #initializing lambda.N at to be consistent with N init
                 s=nimbuild$s,G.true=nimbuild$G.true,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.true3D),
                 y.true=nimbuild$y.true3D,G.latent=nimbuild$G.latent,
                 logit_p0=qlogis(inits$p0),sigma=inits$sigma,lambda=inits$lambda,
                 gammaMat=gammaMat)

#constants for Nimble
constants <- list(M=M,J=J,K=K,K2D=K2D,n.samples=nimbuild$n.samples,n.cat=n.cat,
                xlim=nimbuild$xlim,ylim=nimbuild$ylim,n.levels=n.levels)

#supply data to nimble
Nimdata <- list(y.true=array(NA,dim=c(M,J,K)),
              G.true=matrix(NA,nrow=M,ncol=n.cat),ID=rep(NA,nimbuild$n.samples),
              z=rep(NA,M),X=as.matrix(data$X),capcounts=rep(NA,M))

# set parameters to monitor
parameters <- c('lambda.N','p0','lambda','sigma','N','n','gammaMat')

#can also monitor a different set of parameters with a different thinning rate
parameters2 <- c("ID")
nt <- 5 #thinning rate
nt2 <- 50 #thin more

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
#using config nodes for faster configuration skipping nodes we will assign samplers to below
config.nodes <- c("lambda.N","logit_p0","sigma","lambda","gammaMat")
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt, monitors2=parameters2,thin2=nt2,nodes=config.nodes,
                      useConjugacy = FALSE) 

##Here, we remove the default sampler for y.true
#and replace it with the custom "IDSampler".
#can control ID ups per iteration with IDups argument. But not sure that is always the limiting factor
conf$addSampler(target = paste0("y.true[1:",M,",1:",J,",1:",K,"]"),
                type = 'IDSampler',control = list(M=M,J=J,K=K,this.j=data$this.j,this.k=data$this.k,
                                                  n.samples=nimbuild$n.samples,G.obs=data$G.obs,IDups=1,
                                                  G.obs.seen=G.obs.seen,n.cat=n.cat),
                silent = TRUE)


#replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
#2 options. First one is faster, but does not allow other parameters, say lam0 and sigma, to vary
#as a function of ID covs in G.true
conf$addSampler(target = paste("G.true[1:",M,",1:",n.cat,"]", sep=""),
                type = 'GSampler',
                control = list(M=M,n.cat=n.cat,n.levels=n.levels), silent = TRUE)

z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
#nodes used for update
y.nodes <- Rmodel$expandNodeNames(paste("y.true[1:",M,",1:",J,",1:",K,"]"))
pd.nodes <- Rmodel$expandNodeNames(paste("pd[1:",M,",1:",J,"]"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M,"]"))
calcNodes <- c(N.node,pd.nodes,y.nodes)
conf$addSampler(target = c("N"),
                type = 'zSampler',control = list(z.ups=z.ups,M=M,J=J,
                                                 y.nodes=y.nodes,pd.nodes=pd.nodes,
                                                 N.node=N.node,z.nodes=z.nodes,
                                                 calcNodes=calcNodes),silent = TRUE)


#can set scale and adaptive=FALSE to try not tuning activity centers. no idea what good scale is though. 1/3 sigma seems OK
#activity centers hard to tune, but sSampler only does Metropolis Hastings when z=1, so z=0 updates don't lead to suboptimal
#proposals

for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,scale=1,adaptive=TRUE),silent = TRUE)
}

#use block update for  correlated posteriors. Can use "tries" to control how many times per iteration
conf$addSampler(target = c("logit_p0","sigma","lambda.N"),
                type = 'RW_block',control = list(adaptive=TRUE,tries=1),silent = TRUE)


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
idx <- grep("gamma",colnames(mvSamples)) #can remove gammas from plots if you want to focus on other parameters (can be a lot of them)
plot(mcmc(mvSamples[200:nrow(mvSamples),-idx]))

cor(mcmc(mvSamples[200:nrow(mvSamples),-idx]))
#n is number of individuals captured. True value:
data$n
