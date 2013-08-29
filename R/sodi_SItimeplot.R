#' Plot S and I populations of sodi output over time
#' @import ggplot2 noamtools plyr
#' @export
sodi_SItimepolot <- function(sodi) {
  infpop = ddply(sodi, .(time), function(z) {
             data.frame(time=z$time[1], Infected=factor(c("S","I")), 
                        Count=c(sum(z$alive*(z$infections == 0)),
                                sum(z$alive*(z$infections > 0))))
            })
#  sodi[, Infected:=factor(ifelse(infections == 0 , "S", "I"),
#                          levels=c("S", "I"))]
#  infpop = sodi[, sum(alive), by="time,Infected"]
#  infpop = infpop[order(time, Infected), ]
  plot = ggplot(infpop, aes(x=time, y=Count, fill=Infected, order=Infected)) + 
    geom_area(position="stack", stat="identity") + 
    xlab("Time") + ylab("Population") + theme_nr
  return(plot)
}