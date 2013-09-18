#' Plot S and I populations of sodi output over time
#' @import ggplot2 noamtools plyr
#' @export
sodi_SItimepolot <- function(sodi) {
  infpop = ddply(sodi, .(Time), function(z) {
             data.frame(time=z$Time[1], Infected=factor(c("S","I")), 
                        Count=c(sum(z$Infections == 0),
                                sum(z$Infections > 0)))
            })

  plot = ggplot(infpop, aes(x=time, y=Count, fill=Infected, order=Infected)) + 
    geom_area(position="stack", stat="identity") + 
    xlab("Time") + ylab("Population") + theme_nr
  return(plot)
}