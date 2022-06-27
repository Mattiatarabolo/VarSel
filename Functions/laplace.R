dlaplace <- function(x, mean = 0, rate = 1, log = FALSE){
  par1 <- 1 / (2 * rate)
  par2 <- - abs(x- mean) / rate
  if(log){return(log(par1) + par2)}
  else{return(par1*exp(par2))}
}

rlaplace <- function(n, mean = 0, rate = 1){
  return(-rate*(log(runif(n))-log(runif(n))))
}