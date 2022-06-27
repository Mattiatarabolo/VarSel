# VarSel
R implementation of the environment-dependent birth-death model within a Bayesian framework.

## Functions
Directory [Functions](Functions) contains the main functions that are used in the analyses:
- [halfcauchy](Functions/halfcauchy.R): compute density and random deviates of an half-Cauchy distribution.
- [halfnormal](Functions/halfnormal.R): compute density and random deviates of an half-normal distribution.
- [laplace](Functions/laplace.R): compute density and random deviates of a Laplacian distribution.
- [likelihood_bd_mod_c](Functions/likelihood_bd_mod_c.R): compute the likelihood of a user-defined environment-dependent birth-death process.
- [MLE](Functions/MLE.R): compute the maximum likelihood estimates of the parameters of a user-defined environment-dependent birth-death process.
- [sim_env_bd_mod_c](FUnctions/sim_env_bd_mod_c.R): generate a random phylogenetic tree under a user-defined environment-dependent birth-death process.
- [run_env_bd_MCMC](Functions/run_env_bd_MCMC.R): generate three parallel independent MCMC chains sampled from the posterior of a user-defined environment-dependent birth-death process, using user defined uninformative priors and MH sampling.
- [proposal](Functions/proposal.R): proposal functions for the MH-MCMC.
