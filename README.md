# categorical-SPIM

Nimble MCMC samplers for categorical SPIM from "Spatial capture–recapture for categorically marked populations with an application to genetic capture–recapture"

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.2627

There are two data augmentation approaches, 1) regular data augmentation as used in the original paper and 2) an alternative that allows a Poisson prior on expected abundance (DA2 in file names). The latter allows faster N/z updates.
This is N-prior data augmentation: https://github.com/benaug/SCR-N-Prior-Data-Augmentation

There are 4 observation models, 1) Poisson, 2) Bernoulli, 3) negative binomial, and 4) zero-truncated Poisson hurdle. Currently, 3 is only available with regular data augmentation and 2 and 4 are only available with DA2. The negative binomial requires "better data" to work well, i.e., more partial ID information and/or less home range overlap.

All scripts are set up for >1 ID covariate except there is a Poisson version that allows just 1. Other scripts can be set up for just 1 in the same way.

See test scripts and then the files loaded inside test scrips, e.g., model files, data simulator, data initializer, nimble functions/custom updates.