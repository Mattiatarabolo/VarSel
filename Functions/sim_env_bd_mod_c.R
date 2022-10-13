sim_env_bd_mod <- function(env.func, f.lamb, f.mu, time.stop = 0, max.size = 2000, return.all.extinct = TRUE, prune.extinct = TRUE){
  if (time.stop == 0){stop("Must have stopping time")}
  birthdeath.tree.timevar_simp <- function(f.lamb.env, f.mu.env,time.stop = 0, return.all.extinct = TRUE, 
                                           prune.extinct = TRUE) 
  {
    while (1) {
      nblineages <- c(1)
      times <- c(0)
      b <- f.lamb.env(0)
      d <- f.mu.env(0)
      dt <- rexp(1, (b + d))
      t <- dt
      if (t >= time.stop) {
        t <- time.stop
        alive <- 1
        times <- c(times, t)
        nblineages <- c(nblineages, 1)
        break
      }
      r <- runif(1)
      if (r > b/(b + d)) {
        times <- c(times, dt)
        nblineages <- c(nblineages, 0)
        alive <- rep(FALSE, 1)
      }
      else {
        edge <- rbind(c(1, 2), c(1, 3))
        edge.length <- rep(NA, 2)
        stem.depth <- rep(t, 2)
        alive <- rep(TRUE, 2)
        times <- c(times, dt)
        nblineages <- c(nblineages, sum(alive))
        next.node <- 4
        repeat {
          if (sum(alive) == 0)
            break
          b <- f.lamb.env(t)
          d <- f.mu.env(t)
          dt <- rexp(1, sum(alive) * (b + d))
          t <- t + dt
          if (t >= time.stop || sum(alive) > max.size) {
            t <- time.stop
            times <- c(times, t)
            nblineages <- c(nblineages, sum(alive))
            break
          }
          r <- runif(1)
          if (r <= b/(b + d)) {
            random_lineage <- round(runif(1, min = 1, max = sum(alive)))
            e <- matrix(edge[alive, ], ncol = 2)
            parent <- e[random_lineage, 2]
            alive[alive][random_lineage] <- FALSE
            edge <- rbind(edge, c(parent, next.node), c(parent, next.node + 1))
            next.node <- next.node + 2
            alive <- c(alive, TRUE, TRUE)
            stem.depth <- c(stem.depth, t, t)
            x <- which(edge[, 2] == parent)
            edge.length[x] <- t - stem.depth[x]
            edge.length <- c(edge.length, NA, NA)
            times <- c(times, t)
            nblineages <- c(nblineages, sum(alive))
          }
          else {
            random_lineage <- round(runif(1, min = 1, max = sum(alive)))
            edge.length[alive][random_lineage] <- t - stem.depth[alive][random_lineage]
            alive[alive][random_lineage] <- FALSE
            times <- c(times, t)
            nblineages <- c(nblineages, sum(alive))
          }
        }
      }
      if (return.all.extinct == TRUE | sum(alive) > 0) {
        break
      }
    }
    if (sum(alive) == 0) {
      obj <- NULL
    }
    else if (sum(alive) == 1) {
      obj <- 1
    }
    else {
      edge.length[alive] <- t - stem.depth[alive]
      n <- -1
      for (i in 1:max(edge)) {
        if (any(edge[, 1] == i)) {
          edge[which(edge[, 1] == i), 1] <- n
          edge[which(edge[, 2] == i), 2] <- n
          n <- n - 1
        }
      }
      edge[edge > 0] <- 1:sum(edge > 0)
      tip.label <- 1:sum(edge > 0)
      mode(edge) <- "character"
      mode(tip.label) <- "character"
      obj <- list(edge = edge, edge.length = edge.length,
                  tip.label = tip.label)
      class(obj) <- "phylo"
      obj <- ape::old2new.phylo(obj)
      if (prune.extinct) {
        obj <- geiger::drop.extinct(obj)
      }
    }
    return(list(tree = obj, times = times, nblineages = nblineages))
  }
  f.lamb.env <- function(t) {
    f.lamb(t, env.func(time.stop - t))
  }
  f.mu.env <- function(t) {
    f.mu(t, env.func(time.stop - t))
  }
  res <- birthdeath.tree.timevar_simp(f.lamb.env, f.mu.env, time.stop, return.all.extinct, prune.extinct)
  return(res)
}

sim_env_bd_mod_c <- compiler::cmpfun(sim_env_bd_mod)
rm(sim_env_bd_mod)