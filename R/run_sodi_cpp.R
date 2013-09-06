#' @export
run_sodi2 <- function(parms, init=NULL) {
  init <- if(is.null(inits)) pop <- initiate(parms)
  parms_mod <- parms
  parms_mod$timenames = as.character(times)
  parms_mod <- within(parms_mod, {
         a_ply(c("f", "g", "d", "r", "alpha", "lamda",
                 "beta", "mu", "xi", "omega", "m"),
        1, function(z) {assign(z, c(0, get(z)), envir=sys.frame(-3))})})
  sodi = run_sodi_rcpp(init, parms_mod)
  names(sodi) = parms$times
  sodi = data.table(do.call(rbind,sodi))
  attr(sodi, "parms") = parms
  return(sodi)
}

#' @export
initiate2 <- function(parms) {
  list2env(parms, environment())
  init = data.frame(ID = c(1:n0, rep(0, K-n0)),
                    X = c(runif(n0,min=bbox[1], max=bbox[2]), rep(0, K-n0)),
                    Y = c(runif(n0,min=bbox[3], max=bbox[4]), rep(0, K-n0)),
                    Stage = c(stages0, rep(0, K-n0)),
                    Infections = c(rep(infect0, n0), rep(0, K-n0))
  )
  return(init)

}