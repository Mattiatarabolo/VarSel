#########################                   Logarithmic Psi function            ##############################
.logPsi <- function(t,f.lamb,f.mu,f,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE,dt=0, 
                    r.int.tab = c(NA), r.int.int.tab = c(NA))
{
  if ((cst.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb <- f.lamb(0)
    mu <- f.mu(0)
    r <- lamb-mu
    res <- r*t-2*log(abs(1+(lamb*(exp(r*t)-1)/(r/f))))
    return(res)
  }
  
  ####### exponential dependencies ########
  
  if ((cst.lamb==TRUE) & (expo.mu==TRUE))
    
  {
    lamb0 <- f.lamb(0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(y){lamb0*y-mu0/beta*(exp(beta*y)-1)}
    gvect <- function(y){mapply(r.int,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(y){.Integrate(r.int.0,0,y,stop.on.error=FALSE)}
    res <- r.int(t)-2*log(abs(1+r.int.int(s,t)*f))
    return(res)
  }
  
  if ((expo.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    alpha <- log(f.lamb(1)/lamb0)
    mu0 <- f.mu(0)
    r.int <- function(y){lamb0/alpha*(exp(alpha*y)-1)-mu0*y}
    gvect <- function(y){mapply(r.int,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(y){.Integrate(r.int.0,0,y,stop.on.error=FALSE)}
    res <- r.int(t)-2*log(abs(1+r.int.int(t)*f))
    return(res)
  }
  
  if ((expo.lamb==TRUE) & (expo.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    alpha <- log(f.lamb(1)/lamb0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(y){lamb0/alpha*(exp(alpha*y)-1)-mu0/beta*(exp(beta*y)-1)}
    gvect <- function(y){mapply(r.int,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(y){.Integrate(r.int.0,0,y,stop.on.error=FALSE)}
    res <- r.int(t)-2*log(abs(1+r.int.int(t)*f))
    return(res)
  }
  
  ####### other dependencies ########
  
  else
  {
    if (dt==0)
    {
      # Compute using R integration functions
      r <- function(t){f.lamb(t)-f.mu(t)}
      r.int <- function(y){.Integrate(Vectorize(r),0,y,stop.on.error=FALSE)}
      r.int.0 <- function(y){exp(r.int(y))*f.lamb(y)}
      r.int.int <- function(y){.Integrate(Vectorize(r.int.0),0,y,stop.on.error=FALSE)}
      rst <- r.int(t)
      rist <- r.int.int(t)
      res <- rst-2*log(abs(1+rist*f))
      return(res)
    }
    else
    {
      rst <- r.int.tab[1+floor(t/dt)]
      rist <- r.int.int.tab[1+floor(t/dt)]
      res <- rst-2*log(abs(1+rist*f))
      return(res)
    }
  }
}

.logPsi_c <- compiler::cmpfun(.logPsi)

