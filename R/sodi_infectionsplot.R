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
          geom_bar(stat="identity") + xlim(-1, max(popdist$infections)) +
          scale_y_continuous(limits=c(0, parms$n0), oob=rescale_none) + 
          xlab("Infections") + ylab("Count") + theme_nr
    
  if (length(times) > 1) plot = plot + facet_wrap(~time)
  
  return(plot)

}
