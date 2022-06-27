library(RPANDA)

inputs <- c(commandArgs(trailingOnly=TRUE))
Job <- inputs[1]
which_prior <- "exp"#inputs[2]  #(options "unif" "norm"  "exp")
extinction_rate = "const"#inputs[3]

set.seed(as.numeric(Job))
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


#################                     Phylogenetic tree simulation                   #########################

phylo <- readRDS(paste("/data/biodiv/tarabolo/treesim_twovar/outdata/phylo/phylo_40_", extinction_rate, "_",  Job, ".rds", sep = ""))

# Total time and sampling fraction
tot_time <- max(node.age(phylo)$ages)
tot_time_env = 40
f_total <- 1

# Environmental function definition
load("./Data/InfTemp.R")
load("./Data/co2.R")

env_list <- list("Temperature" = InfTemp, "CO2" = co2)
env_num_tot <- length(env_list)

env_names <- c("Temperature", "CO2")
env_num <- length(env_names)

time_env <- list()
tot_times_env <- list()
env_data <- list()
env_scales <- list()
env_baseline <- list()

for(name in env_names){
  time_env[[name]] <- env_list[[name]][env_list[[name]][,1] <= tot_time_env, 1]
  tot_times_env[[name]] <- max(env_list[[name]][env_list[[name]][,1] <= tot_time_env, 1])
  env_scales[[name]] <-  max((env_list[[name]][env_list[[name]][,1] <= tot_time_env,2]-min(env_list[[name]][env_list[[name]][,1] <= tot_time_env,2])))
  env_baseline[[name]] <- min(env_list[[name]][env_list[[name]][,1] <= tot_time_env,2])
  env_data[[name]] <- (env_list[[name]][env_list[[name]][,1] <= tot_time_env,2] - env_baseline[[name]])/env_scales[[name]]#curves standardization
}

tot_time_env <- Reduce(min, tot_times_env)

if(tot_time_env < tot_time){print("!!!WARNING!!! tot_time_env < tot_phylo_time")}

rm(InfTemp, env_list)

# Smoothing of the env data and environmental function
dof = c("Temperature"=500, "CO2" = 50)

spline_result <- list()
for(j in 1:env_num){
  spline_result[[env_names[j]]] <- smooth.spline(time_env[[env_names[j]]], env_data[[env_names[j]]], df = dof[j])
}

env_func <- function(t) {
  env_pred <- matrix(0, length(t), env_num)
  for (j in 1:env_num){env_pred[,j] = predict(spline_result[[env_names[j]]], t)$y}
  return(env_pred)
}

# In order to perform computation, the env_func is tabulated
# control from lower_bound -10%, upper_bound + 10%
lower_bound_control <- 0.10
upper_bound_control <- 0.10
lower_bound <- as.double(0)
upper_bound <- tot_time
upper_bound_tab <- upper_bound + upper_bound_control*(upper_bound-lower_bound)
lower_bound_tab <- lower_bound - lower_bound_control*(upper_bound-lower_bound)

# Tabulation of the function from lower_bound -10%, to upper_bound + 10%
time_tabulated <- seq(from = lower_bound_tab, to = upper_bound_tab, length.out = 1 + 1e6)
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

# Constant extinction rate, variable speciation rate
f.lamb <- function(t, par){par[1]*exp(env_func_tab(t)%*%par[3:(env_num+2)])}
if (extinction_rate == "const") {
  f.mu <- function(t, par){par[2]}
} else if (extinction_rate == "ratio") {
  f.mu <- function(t, par){par[2]*f.lamb(t,par)}
}


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

par_names <- c("lambda_0", "mu_0", "theta_T", "theta_CO2")

pamhLocalName = paste("/data/biodiv/tarabolo/bayes_bench/outdata/tracer_logfile/tracer_logfile_", extinction_rate, "_", which_prior, "_", Job, sep = "")

sampler <- run_env_bd_MCMC(tree = phylo, f = f_total, f.lamb = f.lamb, f.mu = f.mu, prior = priorDensity, start_gen = parGen, 
                           par_names = par_names, pamhLocalName = pamhLocalName, proposalKernel = "uniform", iteration = 1e4, 
                           thin = 1e2, update = 1e2, adaptation = 1e4, max_iter = 5e5, seed = NULL, nCPU = 3)

#MCMC estimation
chain_list = coda::mcmc.list(lapply(1:3,function(j){coda::mcmc(sampler[[j]]$chain[-(1:10),1:3])}))
par_MCMC = summary(chain_list)

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
  f.lamb.env <- function(t){f.lamb(t, par)}
  f.mu.env <- function(t){f.mu(t, par)}
  ll <- likelihood_bd_mod_c(nbtips, age, tjs, f.lamb.env, f.mu.env, f = 1, dt = 1e-4, cond = "crown")
  return(ll)
}

likelihood <- compiler::cmpfun(likelihood)

par_0_MLE <- rep(0.5, env_num+2)
par_MLE <- MLE(likelihood, par_0_MLE, meth = "Nelder-Mead")

saveRDS(list(MLE = par_MLE, MCMC = par_MCMC, treesize = nbtips), file = paste("/data/biodiv/tarabolo/bayes_bench/outdata/estimates/estimates_", extinction_rate, "_", which_prior, "_", Job, ".rds", sep =  ""))
saveRDS(sampler, file = paste("/data/biodiv/tarabolo/bayes_bench/outdata/sampler/sampler_", extinction_rate, "_", which_prior, "_", Job, ".rds", sep = ""))
