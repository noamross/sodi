#install.packages("~/Dropbox/code/sodi/", type="source", repos=NULL)
library(sodi)
library(manipulate)
library(multicore)
library(spatstat)
library(ggplot2)
library(doMC)
registerDoMC(cores=7)

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
  times = seq(0,300,1), #Times to report.  If a single number, the max time, and all events will be recorded
  lamda_ex = rep(0.00003, 301), #Sequence of external force of infection.  Must be same length as times to report.
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
  
  mg_actions = c(1, rep(0, 300)) #mgmt actions at each time point. Zero for nothing, others are indices of actionlist
  mg_actionlist = c(1)     # 1 for thin_evenly, 2 for thin_spacing
  mg_levels = c(0.5)       # level for each actions
  mg_reprout = c(1)        # allow resprouts after action? 1 for yes, zero for no
  mg_stages = matrix(c(1:4, 0, 0), nrow=1)
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
  times = seq(0,300,1), #Times to report.  If a single number, the max time, and all events will be recorded
  lamda_ex = rep(0.0003, 301), #Sequence of external force of infection.  Must be same length as times to report.
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
  times = seq(0,300,1), #Times to report.  If a single number, the max time, and all events will be recorded
  lamda_ex = rep(0.003, 301), #Sequence of external force of infection.  Must be same length as times to report.
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
  times = seq(0,300,1), #Times to report.  If a single number, the max time, and all events will be recorded
  lamda_ex = rep(0.03, 301), #Sequence of external force of infection.  Must be same length as times to report.
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

sodi = run_sodi(parms=list(vlow=parms1, low=parms2, med=parms3, high=parms4), name="risk", reps=100, progress=TRUE, parallel=TRUE, diagnostics=TRUE)
counts = sodi[Species=="Tanoak", list(Parms=Parms[1], variable=c("T","I"), value=c(sum(Infections >= 0), sum(Infections >=1))), by="Time,Run"]
counts2 = sodi[, list(Parms=Parms[1], Count=sum(Infections >= 0), PctInf=sum(Infections >=1)/sum(Infections >= 0), MeanInf=mean(Infections)), by="Time,Run,Species"]
counts2[, CountL := c(NA,Count[1:(length(Count)-1)]), by="Run,Species"]
counts2[, PctInfL := c(NA,PctInf[1:(length(PctInf)-1)]), by="Run,Species"]
counts2[, TimeL := c(NA, Time[1:(length(Time) -1)]), by ="Run,Species"]
counts2[, MeanInfL := c(NA, MeanInf[1:(length(MeanInf) -1)]), by ="Run,Species"]
counts2 = counts2[!is.na(CountL),]
ggplot(counts2, aes(x=TimeL, y=CountL, xend=Time, yend=Count, col=MeanInfL, group=paste(Run,Species))) + geom_segment(alpha=0.3) + facet_wrap(~Parms) + scale_color_continuous(low="blue", high="red")

sodi <- load_sodi_files("risk-13-10-17-215330")

subsodi <- sodi[Run==20,]
library(manipulate)

sodi_spatialplot(subsodi, 10)

manipulate(sodi_spatialplot(subsodi, TIME), TIME=slider(min=0,max=200,step=1))


ggplot(counts, aes(x=Time, y=value, col=variable, group=paste(Run,variable))) + geom_line(alpha=0.2) + facet_wrap(~Parms)
#calculate mean mortality rate over course of infection
#mortality is alpha times number of infections
mortality = sodi[Species=="Tanoak", list(meanInf = mean(Infections[Infections >= 1]),  
                                         YTD=1/mean(Infections[Infections >= 1]*parms1$alpha[Stage[1]]), 
                                         InfT=sum(Infections >=1)), by="Time,Run,Parms,Stage"]
mortality = mortality[!is.nan(YTD),]
mortality[, drate := 1/YTD]
ggplot(mortality, aes(x=Time, y=meanInf, col=Stage, group=paste(Run,Stage))) + geom_line(alpha=0.5) + facet_wrap(~Parms)
ggplot(subset(mortality, Stage==4), aes(x=InfT, y=YTD, col=as.factor(Stage), group=paste(Run,Stage))) + geom_line(alpha=0.3) + facet_wrap(~Parms)
ggplot(subset(mortality, Stage==4), aes(x=InfT, y=drate, col=as.factor(Stage), group=paste(Run,Stage))) + geom_line(alpha=0.3) + facet_wrap(~Parms)
manipulate(sodi_spatialplot(sodi, TIME), TIME = slider(0, tail(parms$times, 1), step=1))
manipulate(sodi_infectionsplot(sodi, TIME), TIME=slider(0,tail(parms$times,1),1, step=1))
manipulate(sodi_infectionsdensplot(sodi, TIME), TIME=slider(0,tail(parms$times,1),1, step=1))

sodi_SItimepolot(sodi)

#Rprof("out.prof", line.profiling=TRUE)
C_SI = CSI(sodi, progress="time", n.quantiles=6)
#Rprof(NULL)



manipulate(CSI_plot(C_SI, TIME), TIME=slider(0,tail(parms$times,1),0, step=1))
