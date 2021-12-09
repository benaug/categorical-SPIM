# categorical-SPIM

Nimble sampler for categorical SPIM from "Spatial capture–recapture for categorically marked populations with an application to genetic capture–recapture"

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.2627

Currently, only the Poisson observation model is supported without any occasion-level or behavioral response effects. The code is set up for multiple categorical ID covariates, but can be tricked into using just 1 without changing the Nimble code. See test scripts.

The custom nimble functions allow trap effects, though the BUGS code should probably be devectorized across traps when doing that. The categorical partial ID covariates can be used an individual covariates on other parameters, say lam0 or sigma, only if you use "Gsampler2".
