#' Plot sodi output at a fixed point in time.
#' 
#' If there are multiple time points, facet along these
#' @import ggplot2 data.table
#' @export
sodi_spatialplot = function(sodi, times) {
  parms = attr(sodi, "parms")
  max_infections = max(sodi$Infections)  
  n_stages = max(parms$sp_stages)
  
  matched_times = match_times(sodi, times)
  
  plot = ggplot(subset(sodi, Time %in% matched_times),
                aes(x=X, y=Y, size=Stage, col=Infections, shape=Species)) + 
           geom_point() + coord_fixed() + 
           scale_alpha(range=c(0.1,1), limits=c(0,1), breaks=NULL) +
           scale_size(range=c(2,n_stages+1), limits=c(1,n_stages),
                      breaks=1:n_stages) +
           xlim(parms$bbox[1:2]) + ylim(parms$bbox[3:4]) +
           scale_color_gradientn(limits=c(0,max_infections),
                                 colours=c("#408000", "#FFCC66", "#800080"),
                                 values=c(0,1,max_infections)/max_infections)
  
  if (length(times) > 1) plot = plot + facet_wrap(~Time)
    
  return(plot)
}