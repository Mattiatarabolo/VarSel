create_model = function(start_val, par_names, posterior, proposalKernel, tuning=1) {
  model = list()
  model$start_val = start_val
  model$par_names = par_names
  model$posterior = posterior
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
    post = model$posterior(startvalue)
    if (is.infinite(post$LP)) stop("\r","Check if the starting values matches the priors","\r")
    current = post$LP
    mcmc_step[1:(nparam+3)] = c(startvalue, current, post$LL, post$Pr)
    colnames(mcmc_step) = c(model$par_names,"LP","LL","Pr")
    rownames(mcmc_step) = 0
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
    random = ifelse(is.null(param[["random"]]), FALSE, param$random)
    blocks = ifelse(is.null(param[["blocks"]]), 1, param$blocks)
    
    tuning = model$tuning
    
    
    if(length(tuning)!=nparam) tuning = rep(tuning[1],nparam)
    variables = 1:nparam
    
    
    ## --- Prepare the text file (can be used with tracer)
    write.table(data.frame("state"=rownames(mcmc_step),mcmc_step), file=filename, col.names=TRUE, row.names=FALSE, sep = "\t" , 
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
        ## --- randomized indices)
        for(j in variables){
          
          proposalValues = model$proposalfunction2(mcmc_step[1:nparam], tuning = tuning[j], indice=j)
          proposal = proposalValues$proposal
          hastings = proposalValues$hastings
          
          ## ------ Compute the ratio
          
          post_val = model$posterior(proposal)
          newpost = post_val$LP
          probval = newpost - current + hastings
          probab = exp(probval)
          
          # NOTE: there are several waste of times: generation of new proposals, computation of the likelihood for both the tree and priors
          
          ## ------ Evaluate the ratio
          
          if (probval > 0 || probab > runif(1)){
            current = newpost
            mcmc_step[1:nparam] = proposal
            mcmc_step[nparam+1] = newpost
            mcmc_step[nparam+2] = post_val$LL
            mcmc_step[nparam+3] = post_val$Pr
            acc_val[j] = acc_val[j] + 1
          }
        }# End loop over variables
        
        
        
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
      
      ## --- Loop over the variables (Metropolis-within-Gibbs) algorithm
      j = ifelse(random == FALSE, (i-1)%%(nparam)+1, sample(variables,blocks)) 
      # we can put a function outside the loop to avoid the if statement
      
      proposalValues = model$proposalfunction2(mcmc_step[1:nparam], tuning = tuning[j], indice=j)
      proposal = proposalValues$proposal
      hastings = proposalValues$hastings
      
      ## ------ Compute the ratio
      
      post_val = model$posterior(proposal)
      newpost = post_val$LP
      probval = newpost - current + hastings
      probab = exp(probval)
      
      ## ------ Evaluate the ratio
      
      if (probval > 0 || probab > runif(1)){
        current = newpost
        mcmc_step[1:nparam] = proposal
        mcmc_step[nparam+1] = newpost
        mcmc_step[nparam+2] = post_val$LL
        mcmc_step[nparam+3] = post_val$Pr
      }
      
      
      ## ------- Save the chain to a file
      if( i %% thin == 0 ){
        # There is a faster way to save the data?
        write.table(mcmc_step, file=conn, append=TRUE, sep = "\t", col.names=FALSE, row.names=i, quote=FALSE)
      }
      
      if( i %% consoleupdates == 0 & verbose == TRUE){
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
    
    results = list(chain = coda::mcmc(chain), finetune = tuning, acceptance = mcmc_acceptance)
    return(results)
  }


run_env_bd_MCMC = function(tree, f, f.lamb, f.mu, prior = NULL, start_gen, par_names, pamhLocalName, proposalKernel = "bactrian", 
                           iteration = 1e5, thin = 2e3, update = 1e2, adaptation = 1e4, max_iter = 1e6, seed = NULL, nCPU = 3) {
  if(! is.null(seed)) set.seed(seed)
  
  if (!inherits(phylo, "phylo"))
    stop("object \"phylo\" is not of class \"phylo\"")

  npar = length(par_names)

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
    ll <- likelihood_bd_mod_c(nbtips, age, tjs, f.lamb.env, f.mu.env, f, dt = 1e-4, cond = "crown")
    return(ll)
  }
  
  likelihood <- compiler::cmpfun(likelihood)
  
  if (is.null(prior)) prior = function(par)dunif(par, min = c(0, 0, rep(-5, npar)), max = c(5, 5, rep(5,npar)))
    
  post = function(par) {
    if (any(par[1:2] < 0)) {
      return(list(LL = -Inf, LP = -Inf, Pr = -Inf))
    } else {
      LL = likelihood(par)
      Pr = prior(par)
      return(list(LL = LL, LP = LL + Pr, Pr = Pr))
    }
  }
  
  model = create_model(start_val = start_gen, par_names = par_names, posterior = post, proposalKernel = proposalKernel, tuning = 0.1)
  
  ptm = proc.time()
  sampler = parallel::mclapply(1:3, function(j){
    filename=paste(pamhLocalName,"_chain_",j,"_mcmc.log.txt",sep="")
    set.seed(j)
    return(autoMetropolisGibbs(model, iterations = iteration, consoleupdates = 100, thin = thin, autoOptimize = TRUE, 
                                     filename=filename, update = update, adaptation = adaptation, verbose = T))
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
      modelI$start_val = sampler[[j]]$chain[NROW(sampler[[j]]$chain), 1:npar]
      filename=paste(pamhLocalName,"_chain_",j,"_mcmc.log.txt",sep="")
      return(autoMetropolisGibbs(modelI, iterations = iteration, consoleupdates = 100, thin = thin, autoOptimize = F,
                                 filename=filename, update = update, adaptation = adaptation, verbose = T))
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
  
  sampler2 = parallel::mclapply(1:3,function(j){
     set.seed(j)
     modelI = model
     modelI$tuning = sampler[[j]]$finetune
     modelI$start_val = sampler[[j]]$chain[NROW(sampler[[j]]$chain), 1:npar]
     filename=paste(pamhLocalName,"_chain_",j,"_mcmc.log.txt",sep="")
     return(autoMetropolisGibbs(modelI, iterations = iteration, consoleupdates = 100, thin = thin, autoOptimize = F,
                                filename=filename, update = update, adaptation = adaptation, verbose = T))
  }, mc.cores = nCPU, mc.silent = T)

  for(j in 1:3){sampler[[j]]$chain = coda::mcmc(rbind(sampler[[j]]$chain, sampler2[[j]]$chain[-1,]))}
  rep = coda::mcmc.list(lapply(1:3,function(j){coda::mcmc(sampler[[j]]$chain[-(1:10),-c((npar+1):(npar+3))])}))
  gelman = try(coda::gelman.diag(rep))
  if(! inherits(gelman,"try-error")) {gelman = max(gelman$psrf[,2])} else {gelman = 2}
  cat("Final maximum Gelman statistic: ", gelman, "\n")
  ess = try(coda::effectiveSize(rep))
  if(! inherits(ess, "try-error")) {ess = min(ess)} else {ess = 100}
  cat("Final minimum ESS: ", ess, "\n")

  return(sampler)
}
