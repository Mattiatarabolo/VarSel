library(RPANDA)

inputs <- c(commandArgs(trailingOnly=TRUE))
JobId <- inputs[1] # Seed value for the random generation 
which_prior <- inputs[2]  #Options: unif, norm, exp
extinction_rate = inputs[3] #Options: const, ratio

set.seed(as.numeric(JobId))
sessionInfo()

###################                            FUNCTIONS DEFINITION                       ###################
###################                                                                       ###################
#############################################################################################################

source("./Functions/logPhi_c.R")
source("./Functions/likelihood_bd_mod_c.R")
source("./Functions/logPsi_c.R")
source("./Functions/MLE.R")
source("./Functions/halfnormal.R")
source("./Functions/laplace.R")
source("./Functions/run_env_bd_MCMC.R")
source("./Functions/proposal.R")

#############################################################################################################
###################                                                                       ###################
###################                              MAIN  CODE                               ###################
###################                                                                       ###################
#############################################################################################################


#################                   Phylogenetic tree simulation                   #########################

phylo <- readRDS(paste("phylo/phylo_52_", extinction_rate, "_Temperature_",  JobId, ".rds", sep = ""))

# Total time and sampling fraction
tot_time <- max(node.age(phylo)$ages)  #crown age of the phylogeny
f <- 1  #sampling ratio

# Environmental function definition
load("./Data/InfTemp.R")
load("./Data/co2.R")
load("./Data/sealevel.R")
data("silica")
load("./Data/d13c.R")

env_list <- list("Temperature" = InfTemp, "CO2" = co2, "SeaLevel" = sealevel, "d13C" = d13c, "Silica" = silica)
env_list_names = c("Temperature", "CO2", "SeaLevel", "d13C" , "Silica")
env_num_tot <- length(env_list)

env_name <- c("Temperature")
env_num <- 1
tot_time_env = 52  #Used for smoothing of the environmental function

time_env <- list()
env_data <- list()
env_scales <- list()
env_baseline <- list()

for(name in env_name){
  time_env[[name]] <- env_list[[name]][env_list[[name]][,1] <= tot_time_env, 1]
  env_scales[[name]] <-  max((env_list[[name]][env_list[[name]][,1] <= tot_time_env,2]-min(env_list[[name]][env_list[[name]][,1] <= tot_time_env,2])))
  env_baseline[[name]] <- min(env_list[[name]][env_list[[name]][,1] <= tot_time_env,2])
  env_data[[name]] <- (env_list[[name]][env_list[[name]][,1] <= tot_time_env,2] - env_baseline[[name]])/env_scales[[name]]#curves standardization
}

if(tot_time_env < tot_time){print("!!!WARNING!!! tot_time_env < tot_phylo_time")}

rm(InfTemp, env_list, co2, sealevel, silica, d13c)

# Smoothing of the env data and environmental function
dof <- c(500, 30, rep(500, 2), 30) # degrees of freedom for the smoothing of the environmental functions
names(dof) = env_list_names

spline_result <- list()
for(j in 1:env_num){
  spline_result[[env_name[j]]] <- smooth.spline(time_env[[env_name[j]]], env_data[[env_name[j]]], df = dof[j])
}

env_func <- function(t) {
  env_pred <- matrix(0, length(t), env_num)
  for (j in 1:env_num){env_pred[,j] = predict(spline_result[[env_name[j]]], t)$y}
  return(env_pred)
}

# In order to perform computation, the env_func is tabulated
# control from lower_bound -10%, upper_bound + 10%
lower_bound_control <- 0.10
upper_bound_control <- 0.10
lower_bound <- as.double(0)
upper_bound <- tot_time_env
upper_bound_tab <- upper_bound + upper_bound_control*(upper_bound-lower_bound)
lower_bound_tab <- lower_bound - lower_bound_control*(upper_bound-lower_bound)

# Tabulation of the function from lower_bound -10%, to upper_bound + 10%
time_tabulated <- seq(from = lower_bound_tab, to = upper_bound_tab, length.out = 1 + 1e5)
env_tabulated <- env_func(time_tabulated)

# Tabulated function
env_func_tab <- function(t){
  # number of intervals
  n <- NROW(env_tabulated) - 1
  index <- 1 + as.integer( (t - lower_bound_tab) * n / (upper_bound_tab - lower_bound_tab))
  return(matrix(env_tabulated[index,], length(t), env_num))
}

rm(env_data, spline_result, time_env, tot_times_env, dof, j, lower_bound, lower_bound_control, name, time_tabulated,
   upper_bound, upper_bound_control, env_func)

# Choice of the extinction rate
f.lamb <- function(t, par){par[1]*exp(env_func_tab(t)*par[3])}
if (extinction_rate == "const") {
  f.mu <- function(t, par){par[2]}
} else if (extinction_rate == "ratio") {
  f.mu <- function(t, par){par[2]*f.lamb(t, par)}
}

# Priors definition
if (which_prior == "unif") {
  priorDensity <- function(par){
    sum(dunif(par, min = c(rep(0, 2), rep(-5, env_num)), max = rep(5, 2 + env_num), log = T))
  }
  parGen = function(n = 1){
    runif(3, min = c(rep(0, 2), rep(-1, env_num)), max = rep(1, 2 + env_num))
  }
} else if (which_prior == "norm") {
  priorDensity <- function(par){
    sum(dhalfnorm(par[1:2], sigma = rep(2, 2))) + sum(dnorm(par[3:(2 + env_num)], mean = rep(0, env_num), sd = rep(2, env_num)))
  }
  parGen = function(n = 1){
    c(rhalfnorm(2, sigma = rep(1, 2)), rnorm(env_num, mean = rep(0, env_num), sd = rep(1, env_num)))
  }
} else if (which_prior == "exp") {
  priorDensity <- function(par){
    sum(dexp(par[1:2], rate = rep(4, 2))) + sum(dlaplace(par[3:(2 + env_num)], mean = rep(0, env_num), rate = rep(2, env_num)))
  }
  parGen = function(n = 1){
    c(rexp(2, rate = rep(2, 2)), rlaplace(env_num, mean = rep(0, env_num), rate = rep(1, env_num)))
  }
}

priorDensity <- compiler::cmpfun(priorDensity)

par_names <- c("lambda_0", "mu_0", "theta")

# logfile for the MCMC, can be used in tracer
pamhLocalName = paste("tracer_logfile/tracer_logfile_", extinction_rate, "_", which_prior, "_", env_name, "_", JobId, sep = "")

sampler <- run_env_bd_MCMC(tree = phylo, f = f, f.lamb = f.lamb, f.mu = f.mu, prior = priorDensity, start_gen = parGen, 
                           par_names = par_names, pamhLocalName = pamhLocalName, proposalKernel = "uniform", iteration = 1e3, 
                           thin = 1, update = 10, adaptation = 1e2, max_iter = 1e4, seed = NULL, nCPU = 1)

#MCMC estimation
chain_list = coda::mcmc.list(lapply(1:3,function(j){coda::mcmc(sampler[[j]]$chain[-(1:10),1:3])}))
par_MCMC = summary(chain_list)$statistics

#MLE estimation
nbtips <- Ntip(phylo)
from_past <- cbind(phylo$edge,node.age(phylo)$ages)
ages <- rbind(from_past[,2:3],c(nbtips+1,0))
ages <- ages[order(ages[,2]),]
age <- max(ages[,2])
ages[,2] = age- ages[,2]
tjs <- ages[2:(length(ages[,2])-nbtips),2]

rm(from_past, ages)

likelihood <- function(par){
  if (par[1] < 0 || par[2] < 0) return(-Inf)
  f.lamb.env <- function(t){f.lamb(t, par)}
  f.mu.env <- function(t){f.mu(t, par)}
  ll <- likelihood_bd_mod_c(nbtips, age, tjs, f.lamb.env, f.mu.env, f = f, dt = 1e-4, cond = "crown")
  return(ll)
}

likelihood <- compiler::cmpfun(likelihood)

par_0_MLE <- rep(0.5, 3)  #Starting point for the MLE
par_MLE <- MLE(likelihood, par_0_MLE, meth = "Nelder-Mead")

saveRDS(list(MLE = par_MLE, MCMC = par_MCMC, treesize = nbtips), file = paste("bayes_onevar_outdata/estimates_", extinction_rate, "_", which_prior,"_", env_name, "_", Job, ".rds", sep =  ""))
saveRDS(sampler, file = paste("bayes_onevar_outdata/sampler_", extinction_rate, "_", which_prior,"_", env_name, "_", Job, ".rds", sep = ""))