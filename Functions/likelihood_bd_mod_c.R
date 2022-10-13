########################                Logarithmic likelihood             ###################################
likelihood_bd_mod <- function(nbtips, age, tjs, f.lamb, f.mu, f, dt=0, cond = "crown")
{
  
  if(dt != 0)
  {
    X <- seq(0, age + dt, by = dt)
    Nintervals <- length(X)
    r <- function(t){f.lamb(t)-f.mu(t)}
    r.int.tab <- cumsum(r(X)) * dt
    r.int.0 <- function(y){exp(r.int.tab[1 + floor( y/dt )]) * f.lamb(y)}
    r.int.int.tab <- cumsum(r.int.0(X)) * dt
    log_lik_tj <- function(t){
      log(f.lamb(t)) + .logPsi_c(t,f.lamb,f.mu,f,dt=dt, r.int.tab = r.int.tab, r.int.int.tab = r.int.int.tab)}
    log_indLikelihood <- log(f.lamb(age))+2*.logPsi_c(age,f.lamb,f.mu,f,dt=dt, r.int.tab = r.int.tab, r.int.int.tab = r.int.int.tab)
  }
  
  else
  {
    log_lik_tj <- function(t){log(f.lamb(t)) + .logPsi_c(t,f.lamb,f.mu,f,dt=dt)}
    log_indLikelihood <- log(f.lamb(age))+2*.logPsi_c(age,f.lamb,f.mu,f,dt=dt)
  }
  
  vec_lik_tj <- sapply(tjs, log_lik_tj)
  log_indLikelihood <- log_indLikelihood + sum(vec_lik_tj)
  if (!is.finite(log_indLikelihood)){return(-Inf)}
  
  log_data_lik <- log_indLikelihood+nbtips*log(f)
  if (!is.finite(log_data_lik)){return(-Inf)}
  
  # Conditioning
  if (cond==FALSE){log_final_lik <- log_data_lik}
  else{
    if(dt != 0){
      logPhi <- .logPhi_c(age, f.lamb,f.mu,f,dt=dt, r.int.tab = r.int.tab, r.int.int.tab = r.int.int.tab)
    }
    else{
      logPhi <- .logPhi_c(age, f.lamb,f.mu,f, dt=dt)
    }
    if(cond=="stem"){log_final_lik <- log_data_lik - logPhi}
    else if (cond=="crown"){log_final_lik <- log_data_lik - log(f.lamb(age)) - 2*logPhi}
    if (!is.finite(log_final_lik)){return(-Inf)}
  }
  
  return(log_final_lik)
}

likelihood_bd_mod_c <- compiler::cmpfun(likelihood_bd_mod)

rm(likelihood_bd_mod)