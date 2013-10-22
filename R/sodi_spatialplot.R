#' Plot sodi output at a fixed point in time.
#' 
#' If there are multiple time points, facet along these
#' @import ggplot2 data.table plyr
#' @export
sodi_spatialplot = function(sodi, times) {
  parms = attr(sodi, "parms")
  max_infections = max(sodi$Infections)  
  n_stages = length(unique(sodi$Stage))
  matched_times = sodi:::match_times(sodi, times)
  
  
  sodi <- subset(sodi, Time %in% matched_times)
  bbox=c(floor(min(sodi$X)), ceiling(max(sodi$X)), floor(min(sodi$Y)), ceiling(max(sodi$Y)))
  
  sodi$Species <- factor(sodi$Species, levels=c("Tanoak", "Bay", "Redwood"))
  
  sodi <- ddply(sodi, .(Species), function(x) {
    x2 <- x
    if(length(unique(x2$Stage)) == 1) x2$Stage = n_stages
    
    return(x2)
    
  })
  
  plot = ggplot(sodi,
                aes(x=X, y=Y, size=Stage, fill=Infections, col=Species, lwd=Species)) + 
           geom_point(shape=19, aes(size=Stage+3)) +
           geom_point(shape=21, aes(alpha=ifelse(Infections > 0, 1, 0)), oob=rescale_none) + 
           coord_fixed() + 
           scale_alpha(range=c(0.1,1), limits=c(0,1), breaks=NULL) +
           xlim(bbox[1:2]) + ylim(bbox[3:4]) +
           scale_fill_gradientn(limits=c(0,max_infections),
                                 colours=c("#FFFFFF", "#FFCC66", "#990000"),
                                 values=c(0,1,max_infections)/max_infections) +
          scale_size(limits=c(0,n_stages+3),range=c(2,5), guide=FALSE) +
          scale_color_manual(values=c("#408000", "#0000FF","#000000")) + 
          guides(col=guide_legend(override.aes=list(size=5)))
          
  
  if (length(times) > 1) plot = plot + facet_wrap(~Time)
  
  plot <- plot + theme(panel.grid.major.x=element_line(colour="#ECECEC", size=0.5, linetype=1),
                       panel.grid.major.y=element_line(colour="#ECECEC", size=0.5, linetype=1),
                       panel.grid.minor.x=element_blank(),
                       panel.grid.minor.y=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.ticks.y=element_blank(),
                       panel.background=element_blank(),
                       legend.title=element_blank(),
                       legend.key=element_rect(fill="white"),
                       legend.key.size=unit(1.5, "cm"),
                       legend.text=element_text(size=22),
                       axis.title=element_blank(),
                       axis.text=element_blank())
  
  return(plot)
}