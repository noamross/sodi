#' Plot sodi output at a fixed point in time.
#' 
#' If there are multiple time points, facet along these
#' @import ggplot2 data.table
#' @export
sodi_spatialplot = function(sodi, times, dead=FALSE) {
  parms = attr(sodi, "parms")
  max_infections = max(sodi$infections)  
  
  if(!dead) sodi = subset(sodi, alive==1)
  
  matched_times = match_times(sodi, times)
  
  #browser()
  plot = ggplot(subset(sodi, time %in% matched_times),
                aes(x=x, y=y, size=stage, col=infections, alpha=alive)) + 
           geom_point() + coord_fixed() + 
           scale_alpha(range=c(0.1,1), limits=c(0,1), breaks=NULL) +
           scale_size(range=c(2,4), limits=c(1,parms$n_stages),
                      breaks=c(1,parms$n_stages)) +
           xlim(parms$bbox[1:2]) + ylim(parms$bbox[3:4]) +
           scale_color_gradientn(limits=c(0,max_infections),
                                 colours=c("#408000", "#FFCC66", "#800080"),
                                 values=c(0,1,max_infections)/max_infections)
  
  if (length(times) > 1) plot = plot + facet_wrap(~time)
    
#   ggtitle(paste0("Time = ", round(TIME, 1), 
#                  ", Population = ", sum(subst$alive),
#                  ", % Infected = ", round(100*sum(subst$infections > 0 & 
#                                                 subst$alive > 0)/
#                                             sum(subst$alive), 1), "%"))
#  browser()
  return(plot)
}
