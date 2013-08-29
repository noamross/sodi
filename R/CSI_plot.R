#' Plot a heat map of the CSI matrix
#' @import ggplot2 reshape2
#' @export
CSI_plot = function(CSI, times) {
    CSI_times = as.numeric(dimnames(CSI)$time)
    matched_times = CSI_times[approx(x = CSI_times, y = 1:length(CSI_times),
                            xout = times, method = "constant", rule = 2)$y]
    
    d_CSI = melt(CSI[,,CSI_times %in% matched_times])
    lims = c(min(CSI, na.rm=TRUE),max(CSI, na.rm=TRUE))
    lim_rescale = (c(lims[1],0,lims[2]) - lims[1])/(lims[2] - lims[1])
    plot = ggplot(d_CSI, aes(x=Var1, y=Var2, fill=value)) + 
           geom_tile(color="white", lwd=2) +
           coord_fixed() +
           scale_fill_gradientn(limits=lims, colours=c("red", "white", "blue"),
                                values=lim_rescale, space="rgb",
                                na.value="grey80") +
           scale_x_discrete(limits = dimnames(CSI)[[2]] ) +
           scale_y_discrete(limits = dimnames(CSI)[[2]]) + theme_nr +
           theme(panel.grid.major.y=element_blank(),
                 axis.ticks.x=element_blank()) +
           xlab("No. Infections") + ylab("No. Infections")
            
    

    
    if (length(times) > 1) plot = plot + facet_wrap(~time)
    
    return(plot)
    #geom_text - ADD labels on the diagonal of pop count in each bin
  # scale_y_reverse()
}

