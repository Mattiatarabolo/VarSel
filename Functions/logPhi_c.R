#########################                 Logarithmic Phi function              ##############################
.logPhi <- function(t,f.lamb,f.mu,f,cst.lamb=FALSE,cst.mu=FALSE,expo.lamb=FALSE,expo.mu=FALSE,dt=0, 
                    r.int.tab = c(NA), r.int.int.tab = c(NA))
{
  
  if ((cst.lamb==TRUE) & (cst.mu==TRUE))
  {
    lamb <- f.lamb(0)
    mu <- f.mu(0)
    r <- lamb-mu
    res <- log(r)+r*t-log(r/f+lamb*(exp(r*t)-1))
    return(res)
  }
  
  if ((cst.lamb==TRUE) & (expo.mu==TRUE))
  {
    lamb0 <- f.lamb(0)
    mu0 <- f.mu(0)
    beta <- log(f.mu(1)/mu0)
    r.int <- function(y){lamb0*y-mu0/beta*(exp(beta*y)-1)}
    gvect <- function(y){mapply(r.int,y)}
    r.int.0 <- function(y){exp(gvect(y))*f.lamb(y)}
    r.int.int <- function(y){.Integrate(r.int.0,0,y,stop.on.error=FALSE)}
    res <- r.int(t)-log(1/f+r.int.int(t))
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
    res <- r.int(t)-log(1/f+r.int.int(t))
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
    res <- r.int(t)-log(1/f+r.int.int(t))
    return(res)
  }
  
  else
  {
    if (dt==0)
    {
      r <- function(t){f.lamb(t)-f.mu(t)}
      r.int <- function(y){.Integrate(Vectorize(r),0,y,stop.on.error=FALSE)}
      r.int.0 <- function(y){exp(r.int(y))*f.lamb(y)}
      r.int.int <- function(y){.Integrate(Vectorize(r.int.0),0,y,stop.on.error=FALSE)}
      rit <- r.int(t)
      ri0t <- r.int.int(t)
      res <- rit-log(1/f+ri0t)
      return(res)
    }
    else
    {
      rit <- r.int.tab[1+floor(t/dt)]
      ri0t <- r.int.int.tab[1+floor(t/dt)]
      res <- rit-log(1/f+ri0t)
      return(res)
    }
  }
}

.logPhi_c <- compiler::cmpfun(.logPhi)

