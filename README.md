# categorical-SPIM

Nimble sampler for categorical SPIM from "Spatial capture–recapture for categorically marked populations with an application to genetic capture–recapture"

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.2627

Currently, only the Poisson observation model is supported without any occasion or behavioral response effects. The code is set up for multiple categorical ID covariates, but can be tricked into using just 1 without changing the Nimble code. See test scripts.
