require(RPANDA)
library(tictoc)

set.seed(1)
sessionInfo()

inputs <- c(commandArgs(trailingOnly=TRUE))
extinction_rate = "ratio"
#extinction_rate = "ratio"
which_likelihood = inputs#"new"
#which_likelihood = "old"

#############################################################################################################
###################                                                                       ###################
###################                            FUNCTIONS DEFINITION                       ###################
###################                                                                       ###################
#############################################################################################################

source("./Functions/Integrate.R")
source("./Functions/logPhi_c.R")
source("./Functions/likelihood_bd_mod_c.R")
source("./Functions/logPsi_c.R")

#############################################################################################################
###################                                                                       ###################
###################                              MAIN  CODE                               ###################
###################                                                                       ###################
#############################################################################################################
phylo <- readRDS("phylo_40_CO2_1.rds")

if (!inherits(phylo, "phylo"))
  stop("object \"phylo\" is not of class \"phylo\"")

nbtips <- Ntip(phylo)
from_past <- cbind(phylo$edge,node.age(phylo)$ages)
ages <- rbind(from_past[,2:3],c(nbtips+1,0))
ages <- ages[order(ages[,2]),]
age <- max(ages[,2])
ages[,2] = age- ages[,2]
tjs <- ages[2:(length(ages[,2])-nbtips),2]

rm(ages, from_past)

# Total time and sampling fraction
tot_time <- max(node.age(phylo)$ages)
f_total <- 1

# Environmental function definition
load("./Data/InfTemp.R")
load("./Data/co2.R")
load("./Data/sealevel.R")
data("silica")
load("./Data/d13c.R")

env_list <- list("Temperature" = InfTemp, "CO2" = co2, "SeaLevel" = sealevel, "d13C" = d13c, "Silica" = silica)
env_num_tot <- length(env_list)

env_names <- c("Temperature", "CO2" , "SeaLevel", "d13C", "Silica")
env_num <- length(env_names)

time_env <- list()
tot_times_env <- list()
env_data <- list()
env_scales <- list()
env_baseline <- list()

for(name in env_names){
  time_env[[name]] <- env_list[[name]][,1]
  tot_times_env[[name]] <- max(env_list[[name]][,1])
  env_data[[name]] <- (env_list[[name]][,2] - min(env_list[[name]][,2]))/
    max((env_list[[name]][,2]-min(env_list[[name]][,2]))) #curves standardization
  env_scales[[name]] <-  max((env_list[[name]][,2]-min(env_list[[name]][,2])))
  env_baseline[[name]] <- min(env_list[[name]][,2])
}

tot_time_env <- Reduce(min, tot_times_env)

if(tot_time_env < tot_time){print("!!!WARNING!!! tot_time_env < tot_phylo_time")}

rm(InfTemp, co2, sealevel, silica, d13c, env_list, phylo)

# Smoothing of the env data and environmental function
dof = c(rep(50, 4), rep(20, env_num_tot - 4))

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
env_func_tab <- function(t) {
  # number of intervals
  n <- NROW(env_tabulated) - 1
  index <- 1 + as.integer( (t - lower_bound_tab) * n / (upper_bound_tab - lower_bound_tab))
  return(matrix(env_tabulated[index,], length(t), env_num))
}

rm(env_data, spline_result, time_env, tot_times_env, dof, j, lower_bound, lower_bound_control, name, time_tabulated,
   upper_bound, upper_bound_control, env_func)

# Constant extinction rate, variable speciation rate
f.lamb <- function(env, l0, m0, Garray){l0*exp(env%*%Garray)}
if (extinction_rate == "const") {
  f.mu <- function(env, l0, m0, Garray){m0}
} else if (extinction_rate == "ratio") {
  f.mu <- function(env, l0, m0, Garray){m0*f.lamb(env, l0, m0, Garray)}
}


if (which_likelihood == "old") {
  likelihood <- function(l0, m0, Garray, dt){
    f.lamb.env <- function(t){f.lamb(env_func_tab(t), l0, m0, Garray)}
    f.mu.env <- function(t){f.mu(env_func_tab(t), l0, m0, Garray)}
    return(likelihood_bd(phylo, tot_time, f.lamb.env, f.mu.env, f = f_total, cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = FALSE,
                         expo.mu = FALSE, dt = dt, cond = "crown"))
  }
} else if (which_likelihood == "new") {
    likelihood <- function(l0, m0, Garray, dt){
      f.lamb.env <- function(t){f.lamb(env_func_tab(t), l0, m0, Garray)}
      f.mu.env <- function(t){f.mu(env_func_tab(t), l0, m0, Garray)}
      return(likelihood_bd_mod_c(nbtips, age, tjs, f.lamb.env, f.mu.env, f = f_total, dt = dt, cond = "crown"))
  }
}

likelihood_c <- compiler::cmpfun(likelihood)

lls <- matrix(0, nrow = 1000, ncol = 4)
exec_times = matrix(0, nrow = 1000, ncol = 4)

for(i in 1:1000){
  set.seed(i)
  base <- runif(2, min = 0, max = 1)
  Garray <- rnorm(env_num, mean = rep(0, env_num), sd = rep(2, env_num))

  cat("iteration ", i , " dt = 1e-5 ........ ")

  start <- Sys.time()
  lls[i,1] = likelihood_c(base[1], base[2], Garray, dt = 1e-5)
  exec_times[i, 1] = as.numeric(Sys.time() - start)

  cat("DONE\n")
  cat("              dt = 1e-4 ........ ")
  start <- Sys.time()
  lls[i,2] = likelihood_c(base[1], base[2], Garray, dt = 1e-4)
  exec_times[i, 2] = as.numeric(Sys.time() - start)

  cat("DONE\n")
  cat("              dt = 1e-3 ........ ")

  start <- Sys.time()
  lls[i,3] = likelihood_c(base[1], base[2], Garray, dt = 1e-3)
  exec_times[i, 3] = as.numeric(Sys.time() - start)

  cat("DONE\n")
  cat("              dt = 1e-2 ........ ")

  start <- Sys.time()
  lls[i,4] = likelihood_c(base[1], base[2], Garray, dt = 1e-2)
  exec_times[i, 4] = as.numeric(Sys.time() - start)

  cat("DONE\n")
}

saveRDS(list(exectimes = exec_times, likelihoods = lls), paste("/data/biodiv/tarabolo/likelihood_bench/outdata/",which_likelihood,"_bench.rds", sep = ""))
