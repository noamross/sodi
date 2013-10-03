#' @export
#' @import data.table plyr
#' @importFrom plotrix listDepth
run_sodi <- function(parms=NULL, init=NULL, reps=1, progress=FALSE, parallel=FALSE, file=NULL, load=TRUE) {
  if(is.null(file) & load==TRUE) { 
    fname = file.path(tempdir(), "sodi")
  } else if(is.null(file) & load==FALSE) {
    fname = "sodi"
  } else {
    fname = file
  }
  fname = paste0(fname, format(Sys.time(), "%y-%m-%d-%H%M"))
  saveRDS(parms, paste0(fname, ".parmslist.rds"))
  if(listDepth(parms) > 1) {      
    if(is.null(names(parms))) names(parms) <- 1:length(parms)
    reps = rep(reps, length.out=length(parms))
    parms = rep(parms, reps)
    fnames = paste0(fname,".", names(parms), ".", 1:length(parms))
    a_ply(1:length(parms), 1, function(z) {
      run_sodi_single(parms=parms[[z]], init=init, progress=FALSE,
                      filename=fnames[z])
    }, .progress = ifelse(progress, "time", "none"), .parallel=parallel)
                    
    if(load) {
      sodi <- load_sodi_files(fname=fname)
    } else {
      sodi <- data.frame(file=fnames, run=names(parms));
    }
      
  } else {
    run_sodi_single(parms, init=init, progress=progress, filename=name)
    if(load) {
      sodi <- load_sodi_files(fname=fname)
    }
  }
  
  attr(sodi, "parms") = parms
  attr(sodi, "reps") = reps

  
  return(sodi)
}

#' @import stringr
#' @importFrom Hmisc escapeRegex
load_sodi_files <- function(fname, path=".") {
  filenames = list.files(pattern=paste0(escapeRegex(fname), ".*$"), path=path)
  if(length(filenames) <= 2) {
    sodi <- fread(fname)
    parms <- readRDS(paste0(fname,".parmslist.rds"))
    attr(sodi, "parms") <- parms
    return(sodi)
  }
  matches <- str_match(filenames, ".(\\w+).(\\d)$")
  matches <- matches[complete.cases(matches),]
  parms = matches[,2]
  runs = as.numeric(matches[,3])
  sodi <- alply(1:length(filenames), function(z) {
    data <- fread(filenames[z])
    data[, Parms := parms[z]]
    data[, Rep := runs[z]]
  })
  sodi <- rbindlist(sodi)
  attr(sodi, "parms") <- readRDS(paste0(fname, ".parms.rds"))
  return(sodi)
}
    
#' @import data.table plyr
run_sodi_single <- function(parms, init, progress , filename) {
  
  if(is.null(parms)) {
    parms <- sodi::default_parms 
  } else {
    for(name in names(sodi::default_parms)) {
      if(!exists(name, parms))
      parms[[name]] <- sodi::default_parms[[name]]
    }
  }
  
  init <- if(is.null(init)) init <- initiate(parms)

  parms_mod <- parms
  parms_mod$timenames = as.character(parms$times)
  parms_mod <- within(parms_mod, {
         a_ply(c("f", "g", "d", "r", "alpha", "lamda",
                 "beta", "mu", "xi", "omega", "m", "seedm", "max_inf"),
        1, function(z) {assign(z, c(0, get(z)), envir=sys.frame(-3))})})
  sp_stages <- rep(1:length(parms$sp_stages), parms$sp_stages)
  parms_mod$ss <- c(0, rep(which(diff(c(0, sp_stages))==1), parms$sp_stages))
  
  sodi = run_sodi_rcpp(init, parms_mod, progress, filename)
  
  sodi = data.table(rbindlist(sodi))
  
  sodi = sodi[X!=0 & Y!=0 & Stage!=0,]
  sodi[, Species := parms$sp_names[sp_stages[Stage[1]]], by=Stage]
  sodi[, Stage := Stage[1] - parms_mod$ss[Stage[1] + 1] + 1, by=Stage]

}


initiate <- function(parms) {
  list2env(parms, environment())
  init = data.frame(ID = c(1:n0, rep(0, K-n0)),
                    X = c(runif(n0,min=bbox[1], max=bbox[2]), rep(0, K-n0)),
                    Y = c(runif(n0,min=bbox[3], max=bbox[4]), rep(0, K-n0)),
                    Stage = c(stages0, rep(0, K-n0)),
                    Infections = c(rep(1, infect0), rep(0, K-infect0))
                    )
  return(init)
}
  