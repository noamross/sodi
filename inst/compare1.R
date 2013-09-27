#install.packages("~/Dropbox/code/sodi/", type="source", repos=NULL)
library(sodi)
library(manipulate)
library(plyr)
require(ggplot2)

# parms0 <- list(
#   K=100,                 #Carrying capacity
#   bbox = c(0,sqrt(100),0,sqrt(100)),   #Area dimensions
#   n0 = 100,               #Initial population
#   infect0=1,             #Number of infected individuals at start
#   stages0=rep(1,100),     #Distribution of size classes
#   dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
#   seedshadow = 0,          #dispersal kernel for reproduction, if 0, random location
#   sp_names = c("Tanoak"),  #vector of species names
#   sp_stages = c(2),        #vector of the number of size classes for each species
#   m = c(1, 1),           #dispersal parameter
#   seedm = c(1,1),          #dispersal parameter for reproduction
#   times = seq(0,100,1), #Times to report.  If a single number, the max time, and all events will be recorded
#   f = c(0.01, 0.01),     #Fecundity parameter
#   g = c(0.1, 0),         #Growth rates
#   d = c(0.005, 0.005),   #Death rates
#   r = c(0.5, 0.5),       #Resprout probability at death from disease
#   alpha = c(0.005, 0.005),  #Increase in death rate per infection
#   lamda = c(0.3, 0.3),   #Per-stage contact rate.
#   beta = c(1, 1),        #Probability of acquiring disease when contacted
#   mu = c(0.05, 0.05),   #Per-infection recovery rate
#   xi = c(1, 1),     #Fecundity reduction minus one per infection
#   omega = c(1, 1),  #Competitive coefficient
#   max_inf = c(1, 1),  #Maximum number of infections per plant.  Set to 1 for S/I model.
#   beta_meth = 0 #Maximum infection method.  Zero for none, 1 for step function, 2 for decreasing probability
# )
# 
# test0 <- run_sodi(parms0)
# require(microbenchmark)
# test <- microbenchmark(run_sodi(parms0), times=10)

parms1 <- list(
  K=400,                 #Carrying capacity
  bbox = c(0,sqrt(400),0,sqrt(400)),   #Area dimensions
  n0 = 400,               #Initial population
  infect0=1,             #Number of infected individuals at start
  stages0=rep(1,400),     #Distribution of size classes
  dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
  seedshadow = 0,          #dispersal kernel for reproduction, if 0, random location
  sp_names = c("Tanoak"),  #vector of species names
  sp_stages = c(2),        #vector of the number of size classes for each species
  m = c(1, 1),           #dispersal parameter
  seedm = c(1,1),          #dispersal parameter for reproduction
  times = seq(0,250,1), #Times to report.  If a single number, the max time, and all events will be recorded
  f = c(0.01, 0.01),     #Fecundity parameter
  g = c(0.1, 0),         #Growth rates
  d = c(0.005, 0.005),   #Death rates
  r = c(0.5, 0.5),       #Resprout probability at death from disease
  alpha = c(0.005, 0.005),  #Increase in death rate per infection
  lamda = c(0.3, 0.3),   #Per-stage contact rate.
  beta = c(1, 1),        #Probability of acquiring disease when contacted
  mu = c(0.05, 0.05),   #Per-infection recovery rate
  xi = c(1, 1),     #Fecundity reduction minus one per infection
  omega = c(1, 1),  #Competitive coefficient
  max_inf = c(1, 1),  #Maximum number of infections per plant.  Set to 1 for S/I model.
  beta_meth = 0 #Maximum infection method.  Zero for none, 1 for step function, 2 for decreasing probability
)

parms2 <- list(
  K=400,                 #Carrying capacity
  bbox = c(0,sqrt(400),0,sqrt(400)),   #Area dimensions
  n0 = 400,               #Initial population
  infect0=1,             #Number of infected individuals at start
  stages0=rep(1,400),     #Distribution of size classes
  dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
  seedshadow = 0,          #dispersal kernel for reproduction, if 0, random location
  sp_names = c("Tanoak"),  #vector of species names
  sp_stages = c(1),        #vector of the number of size classes for each species
  m = c(1),           #dispersal parameter
  seedm = c(1),          #dispersal parameter for reproduction
  times = seq(0,250,1), #Times to report.  If a single number, the max time, and all events will be recorded
  f = c(0.01),     #Fecundity parameter
  g = c(0),         #Growth rates
  d = c(0.005),   #Death rates
  r = c(0.5),       #Resprout probability at death from disease
  alpha = c(0.005),  #Increase in death rate per infection
  lamda = c(0.3),   #Per-stage contact rate.
  beta = c(1),        #Probability of acquiring disease when contacted
  mu = c(0.05),   #Per-infection recovery rate
  xi = c(1),     #Fecundity reduction minus one per infection
  omega = c(1),  #Competitive coefficient
  max_inf = c(1),  #Maximum number of infections per plant.  Set to 1 for S/I model.
  beta_meth = 0 #Maximum infection method.  Zero for none, 1 for step function, 2 for decreasing probability
)

parms3 <- list(
  K=400,                 #Carrying capacity
  bbox = c(0,sqrt(400),0,sqrt(400)),   #Area dimensions
  n0 = 400,               #Initial population
  infect0=1,             #Number of infected individuals at start
  stages0=rep(1,400),     #Distribution of size classes
  dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
  seedshadow = 3,          #dispersal kernel for reproduction, if 0, random location
  sp_names = c("Tanoak"),  #vector of species names
  sp_stages = c(1),        #vector of the number of size classes for each species
  m = c(1),           #dispersal parameter
  seedm = c(1),          #dispersal parameter for reproduction
  times = seq(0,250,1), #Times to report.  If a single number, the max time, and all events will be recorded
  f = c(0.01),     #Fecundity parameter
  g = c(0),         #Growth rates
  d = c(0.005),   #Death rates
  r = c(0.5),       #Resprout probability at death from disease
  alpha = c(0.005),  #Increase in death rate per infection
  lamda = c(0.3),   #Per-stage contact rate.
  beta = c(1),        #Probability of acquiring disease when contacted
  mu = c(0.05),   #Per-infection recovery rate
  xi = c(1),     #Fecundity reduction minus one per infection
  omega = c(1),  #Competitive coefficient
  max_inf = c(1),  #Maximum number of infections per plant.  Set to 1 for S/I model.
  beta_meth = 0 #Maximum infection method.  Zero for none, 1 for step function, 2 for decreasing probability
)


parms4 <- list(
  K=400,                 #Carrying capacity
  bbox = c(0,sqrt(400),0,sqrt(400)),   #Area dimensions
  n0 = 400,               #Initial population
  infect0=1,             #Number of infected individuals at start
  stages0=rep(1,400),     #Distribution of size classes
  dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
  seedshadow = 3,          #dispersal kernel for reproduction, if 0, random location
  sp_names = c("Tanoak"),  #vector of species names
  sp_stages = c(1),        #vector of the number of size classes for each species
  m = c(1),           #dispersal parameter
  seedm = c(3),          #dispersal parameter for reproduction
  times = seq(0,250,1), #Times to report.  If a single number, the max time, and all events will be recorded
  f = c(0.01),     #Fecundity parameter
  g = c(0),         #Growth rates
  d = c(0.005),   #Death rates
  r = c(0.5),       #Resprout probability at death from disease
  alpha = c(0.005),  #Increase in death rate per infection
  lamda = c(0.3),   #Per-stage contact rate.
  beta = c(1),        #Probability of acquiring disease when contacted
  mu = c(0.05),   #Per-infection recovery rate
  xi = c(1),     #Fecundity reduction minus one per infection
  omega = c(1),  #Competitive coefficient
  max_inf = c(1),  #Maximum number of infections per plant.  Set to 1 for S/I model.
  beta_meth = 0 #Maximum infection method.  Zero for none, 1 for step function, 2 for decreasing probability
)

# require(multicore)
# require(doMC)
# registerDoMC(cores=6)

runs <- run_sodi(parms = list(NoAge=parms1,AgeStruct=parms2,LimitedSeedDisp=parms3,MoreLimitedSeedDisp=parms4), reps=50L, progress=TRUE, parallel=FALSE)
runsum <- runs[, list(Measure=c("Population", "Infected"),
                      Value=c(length(ID), sum(Infections >= 1))),
               by="Time,Parms,Rep"]
runave <- runsum[, list(Value=mean(Value)), by="Time,Parms,Measure"]

plot <- ggplot() + 
  geom_line(data=runave, mapping=aes(x=Time, y=Value, col=Measure), lwd=2) + 
  geom_line(data=runsum, mapping=aes(x=Time, y=Value, col=Measure, group=paste(Rep,Measure)), lwd=0.2, alpha=0.3) + 
  facet_wrap(~Parms)

plot
ggsave("plot.svg", plot)
saveRDS(runs, "test_run_sodi.rds")
# 
# parms4 <- list(
#   K=10000,                 #Carrying capacity
#   bbox = c(0,sqrt(10000),0,sqrt(10000)),   #Area dimensions
#   n0 = 10000,               #Initial population
#   infect0=1,             #Number of infected individuals at start
#   stages0=rep(1,10000),     #Distribution of size classes
#   dispersalfn = 3,        #disease dispersal: 1-exp, 2-fattail, 3-normal,0 for no spatial component
#   seedshadow = 0,          #dispersal kernel for reproduction, if 0, random location
#   sp_names = c("Tanoak"),  #vector of species names
#   sp_stages = c(2),        #vector of the number of size classes for each species
#   m = c(1, 1),           #dispersal parameter
#   seedm = c(1,1),          #dispersal parameter for reproduction
#   times = seq(0,250,1), #Times to report.  If a single number, the max time, and all events will be recorded
#   f = c(0.01, 0.01),     #Fecundity parameter
#   g = c(0.1, 0),         #Growth rates
#   d = c(0.005, 0.005),   #Death rates
#   r = c(0.5, 0.5),       #Resprout probability at death from disease
#   alpha = c(0.005, 0.005),  #Increase in death rate per infection
#   lamda = c(0.3, 0.3),   #Per-stage contact rate.
#   beta = c(1, 1),        #Probability of acquiring disease when contacted
#   mu = c(0.05, 0.05),   #Per-infection recovery rate
#   xi = c(1, 1),     #Fecundity reduction minus one per infection
#   omega = c(1, 1),  #Competitive coefficient
#   max_inf = c(1, 1),  #Maximum number of infections per plant.  Set to 1 for S/I model.
#   beta_meth = 0 #Maximum infection method.  Zero for none, 1 for step function, 2 for decreasing probability
# )
#   
# bigrun <- run_sodi(parms = parms4, progress=TRUE)
# 
