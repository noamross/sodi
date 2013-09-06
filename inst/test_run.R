#install.packages("~/Dropbox/code/sodi/", type="source", repos=NULL)
library(sodi)
library(data.table)
library(spatstat)

parms <- list(
  K=200,                 #Carrying capacity
  bbox = c(0,sqrt(200),0,sqrt(200)),   #Area dimensions
  n0 = 200,               #Initial population
  infect0=1,             #Number of infected individuals at start
  stages0=rep(1,200),     #Distribution of size classes
  dispersal = function(d, a = 1) {
    dnorm(x=d, mean=0, sd=a)
  },
  disp_limit = 10,
  times = seq(0,500,1), #Times to report.  If a single number, the max time, and all events will be recorded
  n_stages = 2,          #Size stages of
  f = c(0.1, 0.1),     #Fecundity parameter
  g = c(0.1, 0),         #Growth rates
  d = c(0.005, 0.005),   #Death rates
  r = c(0.5, 0.5),       #Resprout probability at death from disease
  alpha = c(0.005, 0.005),  #Increase in death rate per infection
  lamda = c(0.3, 0.3),
  beta = c(1, 1),
  mu = c(0.05, 0.05),
  xi = c(1, 1),
  omega = c(1, 1)
)

Rprof("out.prof", line.profiling=TRUE)
sodi = run_sodi(parms, progress=TRUE)
Rprof(NULL)

library(manipulate)
manipulate(sodi_spatialplot(sodi, TIME), TIME = slider(0, tail(parms$times, 1), step=1))

manipulate(sodi_infectionsplot(sodi, TIME), TIME=slider(0,tail(parms$times,1),1, step=1))

manipulate(sodi_infectionsdensplot(sodi, TIME), TIME=slider(0,tail(parms$times,1),1, step=1))

sodi_SItimepolot(sodi)


C_SI = CSI(sodi, progress="time", n.quantiles=6)

manipulate(CSI_plot(C_SI, TIME), TIME=slider(0,tail(parms$times,1),0, step=1))
