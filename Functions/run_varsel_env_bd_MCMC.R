create_model = function(start_val, par_names, posterior, reparametrize, posterior_reparametrize, proposalKernel, tuning=1) {
  model = list()
  model$start_val = start_val
  model$par_names = par_names
  model$posterior = posterior
  model$reparametrize = reparametrize
  model$postrepar = posterior_reparametrize
  model$proposalKernel = proposalKernel
  model$tuning = tuning
  model$proposalfunction2 = function(param, tuning, indice) {
    proposalKernel = model$proposalKernel
    new = proposal_kernel(param[indice], tuning, method = proposalKernel)
    hastings = new$hastings
    param[indice] = new$moves
    
    returnProp = list(
      proposal=param,
      hastings=hastings
    )
    
    return(returnProp)
  }
  return(model)
}

autoMetropolisGibbs = function(model, startvalue, iterations, consoleupdates = 1000, thin = 100, autoOptimize=TRUE, filename,...) {
  
  proposalKernel = model$proposalKernel
  if (missing(startvalue)) {
    startvalue = model$start_val
  }
  
  if (is.function(startvalue)) {
    startvalue = startvalue()
  }
  
  nparam = length(startvalue)
  ## --- Required packages
  # require(data.table)
  ## ---
  mcmc_step = matrix(ncol=nparam+3, nrow=1)
  mcmc_step_transform = matrix(ncol=nparam+3, nrow=1)
  post = model$posterior(startvalue)
  if (is.infinite(post$LP)) stop("\r","Check if the starting values matches the priors","\r")
  current = post$LP
  mcmc_step[1:(nparam+3)] = c(startvalue, current, post$LL, post$Pr)
  mcmc_step_transform[1:(nparam+3)] = c(model$reparametrize(startvalue), current, post$LL, post$Pr)
  colnames(mcmc_step) = c(model$par_names,"LP","LL","Pr")
  colnames(mcmc_step_transform) = c(model$par_names,"LP","LL","Pr")
  rownames(mcmc_step) = 0
  rownames(mcmc_step_transform) = 0
  hastings = 0
  count = 1
  mcmc_acceptance = NULL
  
  if (missing(filename)) {
    filename = paste("mcmc.log.txt")
  }
  
  ## --- Indices for auto-optimization
  kernel = model$proposalKernel
  familyKernel = "uniform"
  if (kernel == "bactrian" | kernel == "Mbactrian" | kernel == "MbactrianLog" | kernel == "MbactrianTriangle" | 
      kernel == "MbactrianLaplace"| kernel == "bactrianTriangle"| kernel == "bactrianLaplace") familyKernel = "bactrian"
  param = list(...)
  Poptimal = ifelse(is.null(param[["acceptance"]]), ifelse(familyKernel == "bactrian", 0.3, 0.44), param$acceptance)
  update = ifelse(is.null(param[["update"]]), 1000, param$update)
  adaptation = ifelse(is.null(param[["adaptation"]]), update*4, param$adaptation)
  verbose = ifelse(is.null(param[["verbose"]]), TRUE, param$verbose)
  blocking = ifelse(is.null(param[["blocking"]]), "no", param$blocking)  # no, random or user
  blocks = ifelse(is.null(param[["blocks"]]), 1, param$blocks)  # if blocking == "no" blocks = 1, if "random" blocks is number of random blocks, if "user" blocks is user-defined list with indices of blocks
  
  tuning = model$tuning
  
  if(length(tuning)!=nparam) tuning = rep(tuning[1],nparam)
  variables = 1:nparam
  
  ## --- Prepare the text file (can be used with tracer)
  write.table(data.frame("state"=rownames(mcmc_step_transform), mcmc_step_transform), file=filename, col.names=TRUE, row.names=FALSE, sep = "\t" , 
              quote=FALSE)
  
  ## --- if autoOptimization==TRUE; we optimize for nrounds during the burnin period
  
  if (autoOptimize == TRUE) {
    nbrounds=adaptation%/%update
    nlbatch=adaptation/nbrounds
    
    # Acceptance matrix
    mcmc_acceptance = matrix(ncol=nparam, nrow=nbrounds)
    
    message("\r","Proposal kernel: ",kernel," with optimal acceptance rate set to: ",Poptimal,"\r")
    message("\r","Optimization of the tuning values for the proposal (",nbrounds," rounds of size: ",nlbatch,") please wait!","\r")
    acc_val = numeric(nparam)
    iterOptim = adaptation
    
    # iteration loop (here it's an adaptation phase we don't records the results!)
    for (i in 1:iterOptim) {
      
      ## --- Loop over the variables (Metropolis-within-Gibbs) algorithm (we can make it random by sampling from a vector of 
      ## --- randomized indices, or from a predefined list)
      if (blocking == "no") {
        indices = variables
      } else if (blocking == "random") {
        indices = sample(variables, blocks)
      } else if (blocking == "user") {
        indices = blocks
      }
      
      for (j in indices) {
        
        proposalValues = model$proposalfunction2(mcmc_step[1:nparam], tuning = tuning[j], indice=j)
        proposal = proposalValues$proposal
        repar_proposal = model$reparametrize(proposal)
        hastings = proposalValues$hastings
        
        ## ------ Compute the ratio
        
        post_val = model$posterior(repar_proposal)
        newpost = post_val$LP
        probval = newpost - current + hastings
        probab = exp(probval)
        post_val_repar = model$posterior(repar_proposal)
        newpost_repar = post_val$LP
        
        # NOTE: there are several waste of times: generation of new proposals, computation of the likelihood for both the tree and priors
        
        ## ------ Evaluate the ratio
        
        if (probval > 0 || probab > runif(1)) {
          current = newpost
          mcmc_step[1:nparam] = proposal
          mcmc_step[nparam+1] = newpost
          mcmc_step[nparam+2] = post_val$LL
          mcmc_step[nparam+3] = post_val$Pr
          mcmc_step_transform[1:nparam] = repar_proposal
          mcmc_step_transform[nparam+1] = newpost_repar
          mcmc_step_transform[nparam+2] = post_val_repar$LL
          mcmc_step_transform[nparam+3] = post_val_repar$Pr
          acc_val[j] = acc_val[j] + 1
        }
      }
      # End loop over variables
      
      ## ------ Update the tuning parameter
      
      if(i %% update == 0){
        Pjump = acc_val/update
        mcmc_acceptance[count,variables] = Pjump
        
        if(any(Pjump < 0.001) == TRUE | any(Pjump > 0.999) == TRUE) {
          indMin = indMax = NULL
          
          if (any(Pjump < 0.001) == TRUE) {
            indMin = which(Pjump < 0.001)
            tuning[indMin] = tuning[indMin]/100
          }
          if (any(Pjump > 0.999) == TRUE) {
            indMax = which(Pjump > 0.999)
            tuning[indMax] = tuning[indMax]*100
          }
          
          indtot = c(indMin, indMax)
          tuning[-indtot] = tuning[-indtot] * (tan((pi/2)*Pjump[-indtot])/tan((pi/2)*Poptimal))
        } else {
          
          # eq. 9 in Yang & Rodriguez 2013 - PNAS
          tuning = tuning * (tan((pi/2)*Pjump)/tan((pi/2)*Poptimal))
        }
        
        # reset the counter
        acc_val[] = 0
        count = count + 1
      }
      
      if( i %% 100 == 0 & verbose==TRUE){
        cat("\r","burnin MCMC in progress",i,"of",iterOptim,"please wait!","\r")
      }
      flush.console()
      
    }
    
    message("\r","Optimization terminated","\r")
    mcmc_acceptance = coda::mcmc(mcmc_acceptance)  # save as a CODA object
  }
  
  
  ## Here if we have optimized the tuning parameters, we continue the mcmc with the burnin values as starting points
  # open a connection to append files
  conn = file(filename, open="a")
  
  for (i in 1:iterations){
    if (blocking == "no") {
      indices = variables
    } else if (blocking == "random") {
      indices = sample(variables, blocks)
    } else if (blocking == "user") {
      indices = blocks
    }
    
    ## --- Loop over the variables (Metropolis-within-Gibbs) algorithm
    for (j in indices) {
      proposalValues = model$proposalfunction2(mcmc_step[1:nparam], tuning = tuning[j], indice=j)
      proposal = proposalValues$proposal
      repar_proposal = model$reparametrize(proposal)
      hastings = proposalValues$hastings
      
      ## ------ Compute the ratio
      
      post_val = model$posterior(repar_proposal)
      newpost = post_val$LP
      probval = newpost - current + hastings
      probab = exp(probval)
      post_val_repar = model$posterior(repar_proposal)
      newpost_repar = post_val$LP
      
      ## ------ Evaluate the ratio
      
      if (probval > 0 || probab > runif(1)){
        current = newpost
        mcmc_step[1:nparam] = proposal  # parameters
        mcmc_step[nparam+1] = newpost
        mcmc_step[nparam+2] = post_val$LL
        mcmc_step[nparam+3] = post_val$Pr
        mcmc_step_transform[1:nparam] = repar_proposal  # parameters
        mcmc_step_transform[nparam+1] = newpost_repar
        mcmc_step_transform[nparam+2] = post_val_repar$LL
        mcmc_step_transform[nparam+3] = post_val_repar$Pr
      }
    }
    
    ## ------- Save the chain to a file
    if( i %% thin == 0 ){
      # There is a faster way to save the data?
      write.table(mcmc_step_transform, file=conn, append=TRUE, sep = "\t", col.names=FALSE, row.names=i, quote=FALSE)
    }
    
    if (i %% consoleupdates == 0 & verbose == TRUE) {
      cat("\r", "MCMC in progress", i, "of", iterations, "please wait!", "\r")
      cat("\n","MCMC in progress", i, "of", iterations, "please wait! \n" )
    } 
    flush.console() 
  }
  
  # close connection
  close(conn)
  message("\r", "Processing the mcmc chain, please wait!", "\r")
  
  
  ## ------- Do the Metropolis-within-Gibbs for separate blocks?
  
  chain = read.table(filename, sep="\t", header=TRUE, row.names=1)
  # To speed up the computations
  #chain = data.table::fread(filename, sep="\t", header=T, skip=1)
  
  results = list(chain = coda::mcmc(chain), finetune = tuning, acceptance = mcmc_acceptance, last_par = mcmc_step[1:nparam])
  return(results)
}


run_varsel_env_bd_MCMC = function(tree, f, f.lamb, f.mu, start_gen, par_names, env_num, pamhLocalName, proposalKernel = "bactrian", 
                           iteration = 1e5, thin = 2e3, update = 1e2, adaptation = 1e4, max_iter = 1e6, seed = NULL, nCPU = 3) 
{
  if (! is.null(seed)) set.seed(seed)
  
  if (!inherits(phylo, "phylo"))
    stop("object \"phylo\" is not of class \"phylo\"")
  
  nbtips <- Ntip(phylo)
  from_past <- cbind(phylo$edge,node.age(phylo)$ages)
  ages <- rbind(from_past[,2:3],c(nbtips+1,0))
  ages <- ages[order(ages[,2]),]
  age <- max(ages[,2])
  ages[,2] = age - ages[,2]
  tjs <- ages[2:(length(ages[,2])-nbtips),2]
  
  npar = length(par_names)
  
  rm(from_past, ages)
  
  likelihood <- function(par){
    if (par[1] < 0 || par[2] < 0) return(-Inf)
    f.lamb.env <- function(t){f.lamb(t, par)}
    f.mu.env <- function(t){f.mu(t, par)}
    ll <- likelihood_bd_mod_c(nbtips, age, tjs, f.lamb.env, f.mu.env, f, dt = 1e-4, cond = "crown")
    return(ll)
  }
  
  likelihood <- compiler::cmpfun(likelihood)
  
  prior = function (par) {
    prior_base = sum(dexp(par[1:2], rate = 3, log = T))
    prior_corr = sum(dnorm(par[3:(2 + env_num)], mean = 0, sd = 1))
    hyperprior_shrink = sum(dunif(par[(3 + env_num):(2 + 2*env_num)], min = 0, max = 1, log = T))
    
    return(prior_base + prior_corr + hyperprior_shrink)
  }
  
  prior = compiler::cmpfun(prior)
  
  prior_repar = function(par) {
    prior_base = sum(dexp(par[1:2], rate = 3, log = T))
    prior_corr = sum(dnorm(par[3:(2 + env_num)], mean = 0, sd = par[(3 + env_num):(2 + 2*env_num)]*par[3 + 2*env_num]))
    hyperprior_shrink = sum(dhalfcauchy(par[(3 + env_num):(3 + 2*env_num)], scale = 1, log = T))
    
    return(prior_base + prior_corr + hyperprior_shrink)
  }
  
  prior_repar = compiler::cmpfun(prior_repar)
  
  post = function(par) {
    if (any(par[c(1:2, (3 + env_num):(3 + 2*env_num))] < 0)) 
    {
      return(list(LL = -Inf, LP = -Inf, Pr = -Inf))
    } else {
      LL = likelihood(par)
      Pr = prior(par)
      return(list(LL = LL, LP = LL + Pr, Pr = Pr))
    }
  }
  
  post_repar = function(par) {
    if (any(par[c(1:2, (3 + env_num):(3 + 2*env_num))] < 0)) 
    {
      return(list(LL = -Inf, LP = -Inf, Pr = -Inf))
    } else {
      LL = likelihood(par)
      Pr = prior_repar(par)
      return(list(LL = LL, LP = LL + Pr, Pr = Pr))
    }
  }
  
  reparametrize = function(par) {
    par_new <- par
    par_new[(3 + env_num):(3 + 2*env_num)] = tan(pi / 2 * par[(3 + env_num):(3 + 2*env_num)])
    par_new[3:(2 + env_num)] = par[3:(2 + env_num)] * par_new[(3 + env_num):(2 + 2*env_num)] * par_new[3 + 2*env_num]
    return(par_new)
  }
  
  model = create_model(start_val = start_gen, par_names = par_names, posterior = post, reparametrize = reparametrize,
                       posterior_reparametrize = post_repar, proposalKernel = proposalKernel, tuning = 0.1)
  
  ptm = proc.time()
  sampler = parallel::mclapply(1:3, function(j){
    filename=paste(pamhLocalName,"_chain_",j,"_mcmc.log.txt",sep="")
    set.seed(j)
    return(autoMetropolisGibbs(model, iterations = iteration, consoleupdates = 100, thin = thin, autoOptimize = TRUE, 
                               filename=filename, update = update, adaptation = adaptation, verbose = T, blocking = "user", 
                               blocks = list(c(1,2, 3 + 2*env_num), c(3: (2+env_num)), c((3 + env_num):(2 + 2*env_num)))))
  }, mc.cores = nCPU, mc.silent = T)
  
  rep = coda::mcmc.list(lapply(1:3,function(j){coda::mcmc(sampler[[j]]$chain[-(1:10),-c((npar+1):(npar+3))])}))
  gelman = try(coda::gelman.diag(rep))
  if(! inherits(gelman,"try-error")) {gelman = max(gelman$psrf[,2])} else {gelman = 2}
  cat("Maximum Gelman statistic: ", gelman, "\n")
  ess = try(coda::effectiveSize(rep))
  if(! inherits(ess, "try-error")) {ess = min(ess)} else {ess = 100}
  cat("Minimum ESS: ", ess, "\n")
  
  tot_iter = iteration
  
  while(((gelman > 1.05) || (ess < 200)) && (tot_iter < max_iter)) {
    sampler2 = parallel::mclapply(1:3,function(j){
      set.seed(j)
      modelI = model
      modelI$tuning = sampler[[j]]$finetune
      modelI$start_val = sampler[[j]]$last_par
      filename=paste(pamhLocalName,"_chain_",j,"_mcmc.log.txt",sep="")
      return(autoMetropolisGibbs(modelI, iterations = iteration, consoleupdates = 100, thin = thin, autoOptimize = F,
                                 filename=filename, update = update, adaptation = adaptation, verbose = T, blocking = "user", 
                                 blocks = list(c(1,2, 3 + 2*env_num), c(3: (2+env_num)), c((3 + env_num):(2 + 2*env_num)))))
    }, mc.cores = nCPU, mc.silent = T)
    
    for(j in 1:3){sampler[[j]]$chain = coda::mcmc(rbind(sampler[[j]]$chain, sampler2[[j]]$chain[-1,]))}
    rep = coda::mcmc.list(lapply(1:3,function(j){coda::mcmc(sampler[[j]]$chain[-(1:10),-c((npar+1):(npar+3))])}))
    gelman = try(coda::gelman.diag(rep))
    if(! inherits(gelman,"try-error")) {gelman = max(gelman$psrf[,2])} else {gelman = 2}
    cat("Maximum Gelman statistic: ", gelman, "\n")
    ess = try(coda::effectiveSize(rep))
    if(! inherits(ess, "try-error")) {ess = min(ess)} else {ess = 100}
    cat("Minimum ESS: ", ess, "\n")
    tot_iter = tot_iter + iteration
  }
  return(sampler)
}