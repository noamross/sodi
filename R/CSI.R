#' Calculate an array of C_SI values.  Takes simulation
#' output of a single time point
#' @import data.table spatstat
#' @export
CSI = function(sodi, times=NULL, progress="none", n.quantiles=6) {
  parms = attr(sodi, "parms")
  if(!is.null(times)) {
    matched_times = match_times(sodi,times)
    sodi = subset(sodi, Time %in% matched_times)
  }
  
  
  #TODO:Make this compatible with a > 1 species, stage scenario.
  parms$dispersalfun <- switch(parms$dispersalfn + 1,
                        function(r) { rep((parms$bbox[2] - parms$bbox[1])*(parms$bbox[4] - parms$bbox[3]), r) },
                        function(r) { expdisp2(r, parms$m[1])[,1] },
                        function(r) { fatdisp2(r, parms$m[1])[,1] },
                        function(r) { normdisp2(r, parms$m[1])[,1] })
                        
  # Assign infection levels into pretty log_binsbins:
  max_infections = max(sodi$Infections)
  breaks = c(0, signif(exp(seq(0, log(max_infections), length.out=n.quantiles)), 
                       digits=1))
  breaks[n.quantiles + 1] = round_any(max_infections,
                                      10^floor(log10(max_infections)), ceiling)
  breaks[breaks > 9] = breaks[breaks > 9] + 1
  
  sodi$mark = cut(sodi$Infections, breaks, include.lowest=TRUE, right=FALSE)
  
  groups = c("0", paste0(breaks[-c(1,n.quantiles+1)], "-", breaks[-c(1,2)] - 1))
  single.bins = which(diff(breaks) == 1)[-1]
  groups[single.bins] = breaks[single.bins]
  levels(sodi$mark) = groups
  
  CSI = aperm(daply(sodi, "Time", function(sodi_step) {
    CSI_step(sodi_step, parms)
    }, .progress=progress), c(2,3,1))
  
  return(CSI)
}

#'Calculate C_SI for a single time step
#'@import spatstat Bolstad
#'@export
CSI_step = function(sodi_step, parms) {
  groupcounts = table(sodi_step$mark)
  groupcounts = groupcounts %o% groupcounts
  sodi_ppp = ppp(x = sodi_step$X, y=sodi_step$Y, 
                 window=owin(parms$bbox[1:2], parms$bbox[3:4]),
                 marks=sodi_step$mark)
  groups = levels(marks(sodi_ppp))
  intensities = intensity(sodi_ppp)
  CSI_s = matrix(NA, length(intensities), length(intensities))
  dimnames(CSI_s) = list(names(intensities),names(intensities))
  oldwarnopt=options()$warn
  options(warn=2)
  for (i in 1:nrow(CSI_s)) {
    for(j in 1:ncol(CSI_s)) {
      if (i >= j) {
        if (groupcounts[i,j] <= 2) next
        PCF = try(pcfcross(sodi_ppp, groups[i], groups[j]), silent=TRUE)
        if("try-error" %in% class(PCF)) next
        disp = parms$dispersalfun(PCF$r)
        fn = disp * (PCF$trans - 1) * PCF$r
        C = suppressMessages(sintegral(PCF$r, fn)$value * 2 * pi)
        CSI_s[i,j] = C
      }
    }
  }
  options(warn=oldwarnopt)
  
  for (i in 1:nrow(CSI_s)) {
    for(j in 1:ncol(CSI_s)) {
      if (i < j) {
        CSI_s[i,j] = CSI_s[j,i]
      }
    }
  }
  return(CSI_s)
}
