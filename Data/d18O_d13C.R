setwd("~/Nextcloud/PersonalFolder/R/Project/mountdir_bioclust/Data")
library(data.table)

d18O.to.Ts <- function(d18O, t){
  d18O.to.Tdo <- function(d18O, t){
    if(0.000 <= t & t < 3.660){return(1 - 4.4 * (( d18O - 3.25) / 3))}
    if(3.600 <= t & t < 34.025){return(5 - 8 * (( d18O - 1.75) / 3))}
    if(34.025 <= t){return(-4*d18O+12)}
  }
  if(0.000 <= t & t < 1.810){return(2*d18O.to.Tdo(d18O, t) + 12.25)}
  if(1.810 <= t & t < 5.330){return(2.5*d18O.to.Tdo(d18O, t) + 12.15)}
  if(5.330 <= t){return(d18O.to.Tdo(d18O, t) + 14.15)}
}

d18O.to.sealevel <- function(d18O){
  if(d18O <= 3.25){return(60 - 40*(d18O - 1.75))}
  if(d18O > 3.25){return(-120*(d18O - 3.25)/1.65)}
}

dT <- function(Ts){return(Ts - 14.15)}

d18O_d13C <- fread("d18O_d13C(2).csv")

time <- d18O_d13C[[1]]

d13C <- d18O_d13C[[2]]

d13c <- data.frame(age = time, d13c = d13C)

save(d13c, file = "d13c.R")

d18O <- d18O_d13C[[3]]

Temp <- mapply(d18O.to.Ts, d18O, time, SIMPLIFY = TRUE)
dTemp <- dT(Temp)

InfTemp <- data.frame(Age = time, Temperature = Temp)

save(InfTemp, file = "InfTemp.R")

seal <- sapply(d18O, d18O.to.sealevel, simplify = TRUE)

sealevel <- data.frame(age = time, sealevel = seal)

save(sealevel, file = "sealevel.R")

co2 <- fread("co2(2).csv")

save(co2, file = "co2.R")
