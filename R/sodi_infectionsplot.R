#' Show the distribution of infection levels
#' @import ggplot2 data.table noamtools
#' @importFrom scales rescale_none
#' @export
sodi_infectionsplot = function(sodi, times) {
  parms = attr(sodi, "parms")
  
  matched_times = match_times(sodi, times)
  
  popdist = sodi[, sum(alive), by="time,infections"]
  plot = ggplot(subset(popdist, time %in% matched_times), 
                aes(x=infections, y=V1)) +
          geom_bar(stat="identity", position="dodge") + xlim(-1, max(popdist$infections)) +
          scale_y_continuous(limits=c(0, parms$n0), oob=rescale_none) + 
          xlab("Infections") + ylab("Count") + theme_nr
    
  if (length(times) > 1) plot = plot + facet_wrap(~time)
  
  return(plot)

}

#' Show the distribution of infection levels
#' @import ggplot2 data.table noamtools
#' @importFrom scales rescale_none
#' @export
sodi_infectionsdensplot = function(sodi, times) {
  parms = attr(sodi, "parms")
  
  matched_times = sodi:::match_times(sodi, times)
  
  popdist = sodi[, sum(alive), by="time,stage,infections"]
  plot = ggplot(subset(popdist, time %in% matched_times), 
                aes(x=infections, fill=as.factor(stage))) +
          geom_density(alpha=0.5, col=0) + xlim(-1, max(popdist$infections)) +
          scale_y_continuous(limits=c(0, 10/max(sodi$infections)), oob=rescale_none) + 
          xlab("Infections") + ylab("Density") + theme_nr
    
  if (length(times) > 1) plot = plot + facet_wrap(~time)
  
  return(plot)
}