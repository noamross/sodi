#install.packages("~/Dropbox/code/sodi/", type="source", repos=NULL)
library(sodi)
library(manipulate)
library(multicore)
library(doMC)
registerDoMC(cores=8)

# Per-plot stem counts
# 1     LIDE         1  5.275229
# 2     LIDE         2 10.018349
# 3     LIDE         3  4.321839
# 4     LIDE         4  3.025316
# 5     LIDE      <NA>  1.000000
# 6     SESE         1  2.357143
# 7     SESE         2  5.542056
# 8     SESE         3  4.435644
# 9     SESE         4  8.147826
# 10    SESE      <NA>  1.000000
# 11    UMCA         1  3.542857
# 12    UMCA         2  4.605634
# 13    UMCA         3  4.647059
# 14    UMCA         4  2.463415
# 15    UMCA      <NA>  1.000000

parms1 <- list(
  K=1168,                 #Carrying capacity
  bbox = c(0,100,0,100),   #Area dimensions
#  n0 = 1000,               #Initial population
  infect0=0,             #Number of infected individuals at start
  randinit=TRUE,
  stages0 = c(0.0106, 0.02, 0.0086, 0.0061, 0.041, 0.0305),
#  stages0=rep(c(1,2,3,4,5,6), c(106,200,86,61,410,305)),     #Distribution of size classes.  n0 length, or, for randiinit, a vector of densities
  dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
  seedshadow = 0,          #dispersal kernel for reproduction, if 0, random location
  sp_names = c("Tanoak", "Bay", "Redwood"),  #vector of species names
  sp_stages = c(4, 1, 1),        #vector of the number of size classes for each species
  m = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),           #dispersal parameter
  seedm = c(1,1, 1, 1, 1, 1),          #dispersal parameter for reproduction
  times = seq(0,100,1), #Times to report.  If a single number, the max time, and all events will be recorded
  lamda_ex = rep(0.00003, 101), #Sequence of external force of infection.  Must be same length as times to report.
  f = c(0, 0.007, 0.2, 0.73, .077, .019),     #Fecundity parameter
  g = c(0.142, 0.2, 0.05, 0, 0, 0),         #Growth rates
  d = c(0.006, 0.003, 0.001, 0.032, .002, .005),   #Death rates
  r = c(0.5, 0.5, 0.5, 0.5, 0, 0),       #Resprout probability at death from disease
  alpha = c(0.005, 0.005, 0.005, .005, 0, 0),  #Increase in death rate per infection
  lamda = c(1, 1, 1, 1, 4, 0),   #Per-stage contact rate.
  beta = c(0.33, 0.32, 0.3, 0.24, 0.3, 0),        #Probability of acquiring disease when contacted
  mu = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05),   #Per-infection recovery rate
  xi = c(1, 1, 1, 1, 1, 1),     #Fecundity reduction minus one per infection
  omega = c(1.93, 0.76, 0.66, 1.55, 1, 1),  #Competitive coefficient
  max_inf =  c(100, 100, 100, 100, 100, 100),  #Maximum number of infections per plant.  Set to 1 for S/I model.
  beta_meth = 2 #Maximum infection method.  Zero for none, 1 for step function, 2 for decreasing probability
)

parms2 <- list(
  K=1168,                 #Carrying capacity
  bbox = c(0,100,0,100),   #Area dimensions
#  n0 = 1000,               #Initial population
  infect0=0,             #Number of infected individuals at start
  randinit=TRUE,
  stages0 = c(0.0106, 0.02, 0.0086, 0.0061, 0.041, 0.0305),
#  stages0=rep(c(1,2,3,4,5,6), c(106,200,86,61,410,305)),     #Distribution of size classes.  n0 length, or, for randiinit, a vector of densities
  dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
  seedshadow = 0,          #dispersal kernel for reproduction, if 0, random location
  sp_names = c("Tanoak", "Bay", "Redwood"),  #vector of species names
  sp_stages = c(4, 1, 1),        #vector of the number of size classes for each species
  m = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),           #dispersal parameter
  seedm = c(1,1, 1, 1, 1, 1),          #dispersal parameter for reproduction
  times = seq(0,100,1), #Times to report.  If a single number, the max time, and all events will be recorded
  lamda_ex = rep(0.0003, 101), #Sequence of external force of infection.  Must be same length as times to report.
  f = c(0, 0.007, 0.2, 0.73, .077, .019),     #Fecundity parameter
  g = c(0.142, 0.2, 0.05, 0, 0, 0),         #Growth rates
  d = c(0.006, 0.003, 0.001, 0.032, .002, .005),   #Death rates
  r = c(0.5, 0.5, 0.5, 0.5, 0, 0),       #Resprout probability at death from disease
  alpha = c(0.005, 0.005, 0.005, .005, 0, 0),  #Increase in death rate per infection
  lamda = c(1, 1, 1, 1, 4, 0),   #Per-stage contact rate.
  beta = c(0.33, 0.32, 0.3, 0.24, 0.3, 0),        #Probability of acquiring disease when contacted
  mu = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05),   #Per-infection recovery rate
  xi = c(1, 1, 1, 1, 1, 1),     #Fecundity reduction minus one per infection
  omega = c(1.93, 0.76, 0.66, 1.55, 1, 1),  #Competitive coefficient
  max_inf =  c(100, 100, 100, 100, 100, 100),  #Maximum number of infections per plant.  Set to 1 for S/I model.
  beta_meth = 2 #Maximum infection method.  Zero for none, 1 for step function, 2 for decreasing probability
)

parms3 <- list(
  K=1168,                 #Carrying capacity
  bbox = c(0,100,0,100),   #Area dimensions
#  n0 = 1000,               #Initial population
  infect0=0,             #Number of infected individuals at start
  randinit=TRUE,
  stages0 = c(0.0106, 0.02, 0.0086, 0.0061, 0.041, 0.0305),
#  stages0=rep(c(1,2,3,4,5,6), c(106,200,86,61,410,305)),     #Distribution of size classes.  n0 length, or, for randiinit, a vector of densities
  dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
  seedshadow = 0,          #dispersal kernel for reproduction, if 0, random location
  sp_names = c("Tanoak", "Bay", "Redwood"),  #vector of species names
  sp_stages = c(4, 1, 1),        #vector of the number of size classes for each species
  m = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),           #dispersal parameter
  seedm = c(1,1, 1, 1, 1, 1),          #dispersal parameter for reproduction
  times = seq(0,100,1), #Times to report.  If a single number, the max time, and all events will be recorded
  lamda_ex = rep(0.003, 101), #Sequence of external force of infection.  Must be same length as times to report.
  f = c(0, 0.007, 0.2, 0.73, .077, .019),     #Fecundity parameter
  g = c(0.142, 0.2, 0.05, 0, 0, 0),         #Growth rates
  d = c(0.006, 0.003, 0.001, 0.032, .002, .005),   #Death rates
  r = c(0.5, 0.5, 0.5, 0.5, 0, 0),       #Resprout probability at death from disease
  alpha = c(0.005, 0.005, 0.005, .005, 0, 0),  #Increase in death rate per infection
  lamda = c(1, 1, 1, 1, 4, 0),   #Per-stage contact rate.
  beta = c(0.33, 0.32, 0.3, 0.24, 0.3, 0),        #Probability of acquiring disease when contacted
  mu = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05),   #Per-infection recovery rate
  xi = c(1, 1, 1, 1, 1, 1),     #Fecundity reduction minus one per infection
  omega = c(1.93, 0.76, 0.66, 1.55, 1, 1),  #Competitive coefficient
  max_inf =  c(100, 100, 100, 100, 100, 100),  #Maximum number of infections per plant.  Set to 1 for S/I model.
  beta_meth = 2 #Maximum infection method.  Zero for none, 1 for step function, 2 for decreasing probability
)

parms4 <- list(
  K=1168,                 #Carrying capacity
  bbox = c(0,100,0,100),   #Area dimensions
#  n0 = 1000,               #Initial population
  infect0=0,             #Number of infected individuals at start
  randinit=TRUE,
  stages0 = c(0.0106, 0.02, 0.0086, 0.0061, 0.041, 0.0305),
#  stages0=rep(c(1,2,3,4,5,6), c(106,200,86,61,410,305)),     #Distribution of size classes.  n0 length, or, for randiinit, a vector of densities
  dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
  seedshadow = 0,          #dispersal kernel for reproduction, if 0, random location
  sp_names = c("Tanoak", "Bay", "Redwood"),  #vector of species names
  sp_stages = c(4, 1, 1),        #vector of the number of size classes for each species
  m = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),           #dispersal parameter
  seedm = c(1,1, 1, 1, 1, 1),          #dispersal parameter for reproduction
  times = seq(0,100,1), #Times to report.  If a single number, the max time, and all events will be recorded
  lamda_ex = rep(0.03, 101), #Sequence of external force of infection.  Must be same length as times to report.
  f = c(0, 0.007, 0.2, 0.73, .077, .019),     #Fecundity parameter
  g = c(0.142, 0.2, 0.05, 0, 0, 0),         #Growth rates
  d = c(0.006, 0.003, 0.001, 0.032, .002, .005),   #Death rates
  r = c(0.5, 0.5, 0.5, 0.5, 0, 0),       #Resprout probability at death from disease
  alpha = c(0.005, 0.005, 0.005, .005, 0, 0),  #Increase in death rate per infection
  lamda = c(1, 1, 1, 1, 4, 0),   #Per-stage contact rate.
  beta = c(0.33, 0.32, 0.3, 0.24, 0.3, 0),        #Probability of acquiring disease when contacted
  mu = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05),   #Per-infection recovery rate
  xi = c(1, 1, 1, 1, 1, 1),     #Fecundity reduction minus one per infection
  omega = c(1.93, 0.76, 0.66, 1.55, 1, 1),  #Competitive coefficient
  max_inf =  c(100, 100, 100, 100, 100, 100),  #Maximum number of infections per plant.  Set to 1 for S/I model.
  beta_meth = 2 #Maximum infection method.  Zero for none, 1 for step function, 2 for decreasing probability
)

sodi = run_sodi(parms=list(vlow=parms1,low=parms2,med=parms3,high=parms4), name="risk", reps=10, progress=TRUE, parallel=TRUE)



manipulate(sodi_spatialplot(sodi, TIME), TIME = slider(0, tail(parms$times, 1), step=1))
manipulate(sodi_infectionsplot(sodi, TIME), TIME=slider(0,tail(parms$times,1),1, step=1))
manipulate(sodi_infectionsdensplot(sodi, TIME), TIME=slider(0,tail(parms$times,1),1, step=1))

sodi_SItimepolot(sodi)

#Rprof("out.prof", line.profiling=TRUE)
C_SI = CSI(sodi, progress="time", n.quantiles=6)
#Rprof(NULL)

manipulate(CSI_plot(C_SI, TIME), TIME=slider(0,tail(parms$times,1),0, step=1))
