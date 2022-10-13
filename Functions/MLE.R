########################                           MLE                     ##################################
MLE <- function(likelihood, par_0, meth = "Nelder-Mead", lower = -Inf)
{
  result_MLE <- suppressWarnings(optim(par_0, likelihood, method = meth, control = list(fnscale = -1), lower = lower))
  return(result_MLE$par)
}
