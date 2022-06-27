# Half-Cauchy or half-t  distribution with 1 degree of freedom
dhalfcauchy <- function(x, scale = 1, log=FALSE){
  if(any(x < 0)) {
    if(log) return(- Inf)
    else return(0)}
  par1 <- 2 / (pi * scale)
  par2 <- 1 / (1 + x^2 / scale^2)
  if(log) {return(log(par1) + log(par2))}
  else return(par1 * par2)
}

rhalfcauchy_unif <- function(n, scale = 1){
  x_tilde <- runif(n)
  x <- scale*tan(pi / 2 * x_tilde)
  return(x)
}

rhalfcauchy_norm <- function(n, scale = 1){
  tau <- rgamma(n, 0.5, 0.5)
  x <- abs(rnorm(n, sd=scale))/sqrt(tau)
  return(x)
}

