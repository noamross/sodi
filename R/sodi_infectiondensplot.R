#' Show the distribution of infection levels
#' @import ggplot2 data.table
#' @importFrom scales rescale_none
#' @export
sodi_infectionsdensplot = function(sodi, times) {
  parms = attr(sodi, "parms")
  
  matched_times = sodi:::match_times(sodi, times)
  
  #popdist = ddply(sodi, .(Time, Stage, Infections), summarize, Count = length(ID))
  popdist = sodi[, length(ID), by="Time,Stage,Infections"]
  setnames(popdist, "V1", "Count")
  
  max_i = max(popdist$Infections)
  max_i = round_any(max_i, 10^(nchar(max_i) - 1), ceiling)
  labels_i = seq(0, max_i, by=10^(nchar(max_i) - 1))
  
  
  plot = ggplot(subset(popdist, Time %in% matched_times), 
                aes(x=Infections, fill=as.factor(Stage))) +
          geom_density(alpha=0.5, col=0) +
          scale_x_continuous(breaks=labels_i, limits=c(0,max_i)) +
          #scale_y_continuous(limits=c(0, 10/max_i), oob=rescale_none) + 
          xlab("Infections") + ylab("Density") # + theme_nr
    
  if (length(times) > 1) plot = plot + facet_wrap(~time)
  
  return(plot)
}