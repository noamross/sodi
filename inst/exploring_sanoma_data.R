require(data.table)
require(ggplot2)
require(plyr)
require(stringr)
require(manipulate)
hs <- read.csv("data/host.stems_sonoma.csv")
hs <- subset(hs, species %in% c("UMCA", "LIDE"))
hs$inf <- ifelse(is.na(hs$symplfct), hs$lidelfct, hs$symplfct)
hs <- subset(hs, !is.na(inf) & inf >= 0)
hs$year <- year(strptime(as.character(hs$date), "%d-%b-%y"))
site = str_match(as.character(hs$plotid), "([[:alpha:]]+)(\\d+)")
hs$site = as.factor(site[,2])
hs$plot = site[,3]
hs <- droplevels(hs)
manipulate(
  ggplot(subset(hs, year==YEAR & species=="UMCA"), aes(x=inf)) + geom_density(fill="red", alpha=0.5) + xlim(0,300) + ylim(0,0.05),
  YEAR=slider(min(hs$year), max(hs$year))
)
  


ggplot(subset(hs, site=="ANN"), aes(x=date, y=symplfct, col=plot, group=tagnumber)) + geom_line(alpha=0.3)
plots <- ddply(hs, .(x,y), function(z) data.frame(x=z$x[1], y=z$y[1], site=z$site[1]))
ggplot(plots, aes(x=x,y=y, col=site)) + geom_point()
hs$year = year(hs$date)
hsa <- hs[, mean(symplfct), by=list(site,year(date))]

