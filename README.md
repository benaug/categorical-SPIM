# categorical-SPIM

Nimble sampler for categorical SPIM from "Spatial capture–recapture for categorically marked populations with an application to genetic capture–recapture"

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.2627

Currently, only the Poisson and Negative Binomial observation models are supported without any occasion-level or behavioral response effects. The code is set up for multiple categorical ID covariates, but can be tricked into using just 1 without changing the Nimble code. See test scripts.

The custom nimble functions allow trap effects, though the BUGS code should probably be devectorized across traps when doing that. The categorical partial ID covariates can be used an individual covariates on other parameters, say lam0 or sigma, only if you use "Gsampler2".

The negative binomial observation model will require "better data" in order to estimate the overdispersion parameter, i.e., more partial ID information and/or less home range overlap. This sampler demonstrates why I've set up the custom ID update the way I have. Using the Poisson observation model, catSPIM can be written in BUGS code without the custom ID update by specifying the *independent* sample-level likelihoods. However, once you switch to any other observation model, you cannot do this. The Metropolis-Hastings update used in the negative binomial sampler will work with any observation model (assuming you switch in the correct likelihood in the custom update).
