library(RPANDA)

sessionInfo()
inputs <- c(commandArgs(trailingOnly=TRUE))
jobId <- "2"#inputs[2]
extinction_rate = "const"#inputs[3]
set.seed(as.numeric(jobId)+1)

#############################################################################################################
###################                                                                       ###################
###################                            FUNCTIONS DEFINITION                       ###################
###################                                                                       ###################
#############################################################################################################
source("./Functions/sim_env_bd_mod_c.R")
source("./Functions/sim_env_bd.R")

#############################################################################################################
###################                                                                       ###################
###################                              MAIN  CODE                               ###################
###################                                                                       ###################
#############################################################################################################


#################                     Phylogenetic tree simulation                   #########################

tot_time <- 40

# Environmental function definition
load("./Data/InfTemp.R")
load("./Data/co2.R")
load("./Data/sealevel.R")
data("silica")
load("./Data/d13c.R")

env_list <- list("Temperature" = InfTemp, "CO2" = co2, "SeaLevel" = sealevel, "d13C" = d13c, "Silica" = silica)
env_list_names = c("Temperature", "CO2", "SeaLevel", "d13C" , "Silica")
env_num_tot <- length(env_list)

env_name <- "Temperature"#inputs[1]
env_num <- 1

time_env <- list()
tot_times_env <- list()
env_data <- list()
env_scales <- list()
env_baseline <- list()

for(name in env_name){
  time_env[[name]] <- env_list[[name]][env_list[[name]][,1] <= tot_time, 1]
  tot_times_env[[name]] <- max(env_list[[name]][env_list[[name]][,1] <= tot_time, 1])
  env_scales[[name]] <-  max((env_list[[name]][env_list[[name]][,1] <= tot_time,2]-min(env_list[[name]][env_list[[name]][,1] <= tot_time,2])))
  env_baseline[[name]] <- min(env_list[[name]][env_list[[name]][,1] <= tot_time,2])
  env_data[[name]] <- (env_list[[name]][env_list[[name]][,1] <= tot_time,2] - env_baseline[[name]])/env_scales[[name]]#curves standardization
}

tot_time_env <- Reduce(min, tot_times_env)

if(tot_time_env < tot_time){print("!!!WARNING!!! tot_time_env < tot_phylo_time")}

rm(InfTemp, co2, sealevel, silica, d13c)

# Smoothing of the env data and environmental function
dof <- c(500, 100, rep(500, 2), rep(50, env_num_tot - 4))
names(dof) = env_list_names

spline_result <- list()
for(name in env_name){
  spline_result[[name]] <- smooth.spline(time_env[[name]], env_data[[name]], df = dof[[name]])
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

rm(env_data, env_list, spline_result, time_env, tot_times_env, dof, lower_bound, lower_bound_control, name, time_tabulated,
   upper_bound, upper_bound_control, env_func)

par_names <- c("lambda_0", "mu_0",
               "G_T", "G_CO2", "G_sea", "G_Si", "G_d13C")

if (extinction_rate == "const") {
    par <- c(0.04, 0.02, 0)
  } else if (extinction_rate == "ratio") {
    par <- c(0.02, 0.02, 1.1)
  }


# Constant extinction rate, variable speciation rate
f.lamb <- function(t, env){par[1]*exp(env*par[3])}
if (extinction_rate == "const") {
  f.mu <- function(t, env){par[2]}
} else if (extinction_rate == "ratio") {
  f.mu <- function(t, env){par[2]*f.lamb(t, env)}
}

set.seed(3)
phyloext <- NULL
while(class(phyloext) != "phylo"){
  phyloext <- sim_env_bd_mod_c(env.func = env_func_tab, f.lamb = f.lamb, f.mu = f.mu, time.stop = tot_time, return.all.extinct = TRUE,
                            prune.extinct = F)$tree
}
phyloext$edge.length[5] = 2
phyloext$edge.length[7] = 4
phyloext$edge.length[8] = phyloext$edge.length[8] -2
phyloext$edge.length[4] = phyloext$edge.length[4] - 2
phyloext$edge.length[9] = phyloext$edge.length[9] + 2
phyloext$edge.length[10] = phyloext$edge.length[10] + 2
plot(phyloext, type = "phylogram", node.pos = NULL, show.tip.label = F, show.node.label = F, 
     edge.color = c(rep("red", 6), "black", rep("red", 2), "black"), edge.width= c(rep(3, 6), 1, rep(3, 2), 1), plot = TRUE)
axisPhylo()
mtext("Time (Myrs)", side = 1, line = 3.3, at = 10, cex = 2.5)


set.seed(3)
phylo <- NULL
while(class(phylo) != "phylo"){
  phylo <- sim_env_bd_mod_c(env.func = env_func_tab, f.lamb = f.lamb, f.mu = f.mu, time.stop = tot_time, return.all.extinct = TRUE,
                            prune.extinct = T)$tree
}
phylo$edge = phylo$edge[-7,]
phylo$edge.length = phylo$edge.length[-7]
#phylo$Nnode = phylo$Nnode - 1
par(cex.axis=2.5, cex.lab=2.5, cex.main=2.5, cex.sub=2.5)
# phylo$tip.label  =as.character(1:4)
plot(phylo, type = "phylogram", node.pos = NULL, show.tip.label = F, show.node.label = F,
     edge.color = NULL, edge.width= 3, edge.lty = NULL, node.color = "red", node.width = 3, tip.color = "red", 
     plot = TRUE)
axisPhylo()
mtext("Time (Myrs)", side = 1, line = 3.3, at = 10, cex = 2.5)

#saveRDS(phylo, file = paste("/data/biodiv/tarabolo/treesim_onevar/outdata/", extinction_rate, "/phylo_40_", env_name, "_", jobId, ".rds", sep = ""))
#saveRDS(par, file = paste("/data/biodiv/tarabolo/treesim_onevar/outdata/par/par_40_", env_name, "_", jobId, ".rds", sep = ""))