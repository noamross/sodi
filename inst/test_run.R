#install.packages("~/Dropbox/code/sodi/", type="source", repos=NULL)
library(sodi)
library(data.table)
library(spatstat)
library(plyr)

parms <- list(
  K=200,                 #Carrying capacity
  bbox = c(0,sqrt(200),0,sqrt(200)),   #Area dimensions
  n0 = 200,               #Initial population
  infect0=2,             #Number of infected individuals at start
  stages0=rep(1,200),     #Distribution of size classes
  dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
  seedshadow = 0,          #dispersal kernel for reproduction, if 0, random location
  sp_names = c("Tanoak"),  #vector of species names
  sp_stages = c(2),        #vector of the number of size classes for each species
  m = c(1, 1),           #dispersal parameter
  seedm = c(1,1),          #dispersal parameter for reproduction
  times = seq(0,100,1), #Times to report.  If a single number, the max time, and all events will be recorded
  n_stages = 2,          #Size stages of
  f = c(0.1, 0.1),     #Fecundity parameter
  g = c(0.1, 0),         #Growth rates
  d = c(0.005, 0.005),   #Death rates
  r = c(0.5, 0.5),       #Resprout probability at death from disease
  alpha = c(0.005, 0.005),  #Increase in death rate per infection
  lamda = c(0.3, 0.3),   #Per-stage contact rate.
  beta = c(1, 1),        #Probability of acquiring disease when contacted
  mu = c(0.05, 0.05),   #Per-infection recovery rate
  xi = c(1, 1),     #Fecundity reduction minus one per infection
  omega = c(1, 1),  #Competitive coefficient
  max_inf = 999999  #Maximum number of infections per plant.  Set to 1 for S/I model.
)

sodi = run_sodi(parms, progress=TRUE)


library(manipulate)

manipulate(sodi_spatialplot(sodi, TIME), TIME = slider(0, tail(parms$times, 1), step=1))
manipulate(sodi_infectionsplot(sodi, TIME), TIME=slider(0,tail(parms$times,1),1, step=1))
manipulate(sodi_infectionsdensplot(sodi, TIME), TIME=slider(0,tail(parms$times,1),1, step=1))

sodi_SItimepolot(sodi)

#Rprof("out.prof", line.profiling=TRUE)
C_SI = CSI(sodi, progress="time", n.quantiles=6)
#Rprof(NULL)

manipulate(CSI_plot(C_SI, TIME), TIME=slider(0,tail(parms$times,1),0, step=1))
