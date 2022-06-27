dhalfnorm <- function(x, sigma = 1, log = FALSE){
  if(any(x < 0)) {
    if(log) return(- Inf)
    else return(0)}
  par1 <- sqrt(2 / pi) / sigma
  par2 <- - x^2 / (2 * sigma)
  if(log){return(log(par1) + par2)}
  else{return(par1*exp(par2))}
}

rhalfnorm <- function(n, sigma = 1){
  return(abs(rnorm(n, mean = 0, sd = sigma)))
}