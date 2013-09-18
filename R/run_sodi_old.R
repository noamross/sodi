#' Run the individual-based SOD model
#' @import data.table
run_sodi_old <- function(parms, progress=FALSE, inits=NULL) {
  if(is.null(inits)) pop <- initiate_old(parms)
  list2env(parms, environment())
  time_max = tail(times, 1) 
  count = n0
  time = 0
  step = 0
  S = pop[,"stage"]
  I = pop[,"infections"]
  E = 1 - sum(pop[, "alive"] * omega[S])/K
  irates = E*f[S]*xi[S]^I + d[S] + g[S] + pop[,"force"] + I * mu[S] + I * alpha[S]
  sodi = list()
  j <- rep(NA, ceiling(sum(irates) * time_max))
  actions <- j
  sodi[[1]] = pop
  
  if(progress) {
    pb = txtProgressBar(0, time_max, 0, style=3)
  }
  while(time < time_max) {
    step = step + 1
    last_time = time
    time = time + rexp(1, sum(irates))
    if(length(times) > 1) {
      if(sum(times > last_time & times < time) > 0) {
        times[times > last_time & times < time][1]
        sodi[[step + 1]] = pop
        sodi[[step + 1]][, "time"] = times[times > last_time & times < time][1]
      }
    } else {
      sodi[[step + 1]] = pop
    }
    
    j = sample(nrow(pop),1, prob=irates)
    s = pop[j,"stage"]
    i = pop[j,"infections"]
    al = pop[j, "alive"]
    pop[, "time"] = time
    possibilities = c(
      reproduce = unname(al*E*f[s]*xi[s]^i),
      die = unname(al*d[s]),
      grow = unname(al*g[s]),
      sicken = unname(al*pop[j,"force"]),
      recover = unname(al*pop[j, "infections"] * mu[s]),
      succumb = unname(al*pop[j, "infections"] * alpha[s] * (1 - r[s])),
      resprout = unname(al*pop[j, "infections"] * alpha[s] * r[s])
    )
    action = sample(names(possibilities), 1, prob=possibilities)
    pop[,"force"] = pmax(1e-99, pop[,"force"] + switch(action,
                                                       reproduce=0,
                                                       die = -force_from(pop,j, parms),
                                                       grow = -force_from(pop,j, parms) + force_from_grow(pop,j, parms),
                                                       sicken = force_from_one(pop, j, parms),
                                                       recover = -force_from_one(pop, j, parms),
                                                       succumb = -force_from(pop, j, parms),
                                                       resprout = -force_from(pop, j, parms)                                    
    ))
    pop = do.call(action, list(pop, j, parms, count, step, time))
    
    S = pop[,"stage"]
    I = pop[,"infections"]
    E = max(0, 1 - sum(pop[, "alive"] * omega[S])/K)
    irates = pop[,"alive"]*(E*f[S]*xi[S]^I + d[S] + g[S] + 
                              pop[,"force"] + I * mu[S] + I * alpha[S])
    if(sum(irates) <= 0) time = time_max 
    if (progress) {
      setTxtProgressBar(pb, time)
    }
  }
  sodi = data.table(do.call(rbind,sodi))
  attr(sodi, "parms") = parms
  return(sodi)
}


#'Set initial conditions and vectors
#' @export
initiate_old <- function(parms) {
  list2env(parms, environment())
  pop =  matrix(NA, nrow=n0, ncol=9)
  colnames(pop) = c("step", "time", "ID","x","y","stage","infections", "alive", "force")
  pop[,"step"] = 0
  pop[,"time"] = 0
  pop[,"ID"] = 1:n0
  pop[,"x"] = runif(nrow(pop),min=bbox[1], max=bbox[2])
  pop[,"y"] = runif(nrow(pop),min=bbox[1], max=bbox[2])
  pop[,"stage"] = stages0
  pop[,"alive"] = 1
  pop[,"infections"] = 0
  pop[1:infect0,"infections"] = 1
  pop[,"force"] = force_on(pop, 1:nrow(pop), parms)
  return(pop)
}

#' @import plyr
force_on <- function(pop, j, parms) {
  aaply(j, 1, function(z) {
    parms$beta[pop[z,"stage"]] * sum(parms$dispersal(distance(pop[z, "x"]-pop[, "x"], pop[z, "y"]-pop[, "y"])) * pop[,"infections"] * parms$lamda[pop[,"stage"]] * pop[, "alive"])
  })
}

force_from <- function(pop, j, parms) {
  parms$beta[pop[,"stage"]] * parms$dispersal(distance(pop[, "x"] - pop[j, "x"], pop[, "y"] - pop[j, "y"])) * pop[j, "infections"] * parms$lamda[pop[j, "stage"]] * pop[j, "alive"]
}

force_from_one <- function(pop, j, parms) {
  parms$beta[pop[,"stage"]] * parms$dispersal(distance(pop[, "x"] - pop[j, "x"], pop[, "y"] - pop[j, "y"])) * parms$lamda[pop[j, "stage"]] * pop[j, "alive"]
}

force_from_grow <- function(pop, j, parms) {
  parms$beta[pop[,"stage"]] * parms$dispersal(distance(pop[, "x"] - pop[j, "x"], pop[, "y"] - pop[j, "y"])) * pop[j, "infections"]  * parms$lamda[pop[j, "stage"] + 1] * pop[j, "alive"]
}

reproduce = function(pop, j, parms, count, step, time) {
  count = count + 1
  pop = rbind(pop, c(step, time, count, runif(1, min=parms$bbox[1], max=parms$bbox[2]),
                     runif(1, min=parms$bbox[3], max=parms$bbox[4]),
                     1, 0, 1, 0))
  pop[count, "force"] = force_on(pop, count, parms)
  assign("count", count, parent.frame())
  return(pop)
}

die = function(pop, j, parms, ...) {
  pop[j, "alive"] = 0
  return(pop)
}

grow = function(pop, j,...) {
  pop[j, "stage"] = pop[j, "stage"] + 1
  return(pop)
}

sicken = function(pop, j,...) {
  pop[j, "infections"] = pop[j, "infections"] + 1
  return(pop)
}

recover = function(pop, j,...) {
  pop[j, "infections"] = pop[j, "infections"] - 1
  return(pop)
}

succumb = function(pop, j, parms, ...) {
  return(die(pop, j, parms))
}

resprout = function(pop, j, ...) {
  pop[j, "stage"] = 1
  pop[j, "infections"] = 0
  return(pop)
}


distance <- function(x,y) {
  sqrt(x^2 + y^2)
}