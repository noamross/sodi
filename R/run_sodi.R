#' @export
#' @import data.table plyr
run_sodi <- function(parms, init=NULL, progress=FALSE) {
  init <- if(is.null(init)) init <- initiate(parms)
  
  parms_mod <- parms
  parms_mod$timenames = as.character(parms$times)
  parms_mod <- within(parms_mod, {
         a_ply(c("f", "g", "d", "r", "alpha", "lamda",
                 "beta", "mu", "xi", "omega", "m", "seedm"),
        1, function(z) {assign(z, c(0, get(z)), envir=sys.frame(-3))})})
  sp_stages <- rep(1:length(parms$sp_stages), parms$sp_stages)
  parms_mod$ss <- c(0, rep(which(diff(c(0, sp_stages))==1), parms$sp_stages))
  
  sodi = run_sodi_rcpp(init, parms_mod, progress)
  
  sodi = data.table(do.call(rbind,sodi))
  
  sodi = sodi[X!=0 & Y!=0 & Stage!=0,]
  sodi[, Species := parms$sp_names[sp_stages[Stage[1]]], by=Stage]
  sodi[, Stage := Stage[1] - parms_mod$ss[Stage[1] + 1] + 1, by=Stage]
  attr(sodi, "parms") = parms
  
  return(sodi)
}

initiate <- function(parms) {
  list2env(parms, environment())
  init = data.frame(ID = c(1:n0, rep(0, K-n0)),
                    X = c(runif(n0,min=bbox[1], max=bbox[2]), rep(0, K-n0)),
                    Y = c(runif(n0,min=bbox[3], max=bbox[4]), rep(0, K-n0)),
                    Stage = c(stages0, rep(0, K-n0)),
                    Infections = c(rep(1, infect0), rep(0, K-infect0))
                    )
  return(init)
}