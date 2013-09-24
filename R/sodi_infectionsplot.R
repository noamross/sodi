#' Show the distribution of infection levels
#' @import ggplot2 data.table plyr
#' @importFrom scales rescale_none
#' @export
sodi_infectionsplot = function(sodi, times) {
  parms = attr(sodi, "parms")
  
  matched_times = sodi:::match_times(sodi, times)
  
  
  popdist = sodi[, length(ID), by="Time,Infections"]
  setnames(popdist, "V1", "Count")
  #popdist = ddply(sodi, .(Time, Stage, Infections), summarize, Count = length(ID))
  
  max_i = max(popdist$Infections)
  max_i = round_any(max_i, 10^(nchar(max_i) - 1), ceiling)
  labels_i = seq(0, max_i, by=10^(nchar(max_i) - 1))
  
  plot = ggplot(subset(popdist, Time %in% matched_times), 
                aes(x=Infections, y=Count)) +
          geom_bar(stat="identity", position="dodge", width=1) +
          scale_x_continuous(breaks=labels_i, limits=c(-0.5,max_i + 0.5)) +
          scale_y_continuous(limits=c(0, parms$n0), oob=rescale_none) + 
          xlab("Infections") + ylab("Count") #+ theme_nr
    
  if (length(times) > 1) plot = plot + facet_wrap(~time)
  
  return(plot)

}

