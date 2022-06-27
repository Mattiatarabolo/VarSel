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

## Data
Directory [Data](Data) constains the environmental variables used in the analyses.

## Phylo
Directory [phylo](phylo) constain the phylogenetic trees used in the analyses.

## [Treesim_onevar](Treesim_onevar.R)
Example of how simulating a phylogenetic tree using only one environmental dependency.

```
Rscript Treesim_onevar.R [env_var] [seed] [extinction_rate]
```
Options for `env_var`: `Temperature`, `CO2`, `SeaLevel`, `d13C`, `Silica`.
Options for `extinction_rate`: `const`, `ratio`.

## [Treesim_onevar](Treesim_onevar.R)
Example of how simulating a phylogenetic tree using two environmental dependencies.

```
Rscript Treesim_onevar.R [env_var_1] [env_var_2] [seed] [extinction_rate]
```
Options for `env_var_1` and `env_var_2`: `Temperature`, `CO2`, `SeaLevel`, `d13C`, `Silica`.
Options for `extinction_rate`: `const`, `ratio`.

## [Bayes_onevar](Bayes_onevar.R)
Example of how performing a MH-MCMC sampling of the posterior with only one environmental dependency.

```
Rscript Bayes_onevar.R [seed] [prior] [extinction_rate]
```
Options for `prior`: `unif`, `norm`, `exp`.
Options for `extinction_rate`: `const`, `ratio`.

## [Bayes_twovar](Bayes_twovar.R)
Example of how performing a MH-MCMC sampling of the posterior with two environmental dependencies.

```
Rscript Bayes_onevar.R [seed] [prior] [extinction_rate]
```
Options for `prior`: `unif`, `norm`, `exp`.
Options for `extinction_rate`: `const`, `ratio`.