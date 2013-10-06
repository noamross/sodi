#' @export
#' @import data.table plyr
#' @importFrom plotrix listDepth
run_sodi <- function(parms=NULL, init=NULL, reps=1, progress=FALSE, parallel=FALSE, name=NULL, load=TRUE, append_date=TRUE) {
  START = Sys.time()
  wd = getwd()
  if (is.null(name)) {
    fname="sodi"
  } else {
    fname=name
  }
  if (append_date==TRUE) fname = paste0(fname, format(Sys.time(), "-%y-%m-%d-%H%M%S"))
  if (is.null(name) & load==TRUE) {
    dir = file.path(tempdir(), fname) 
  } else {
    dir = file.path(wd, fname)
  }
  dir.create(dir)
  setwd(dir)
  if(listDepth(parms) > 1 || reps > 1) {
    if(listDepth(parms) == 1) parms=list(parms)
    if(is.null(names(parms))) names(parms) <- 1:length(parms)
    reps = rep(reps, length.out=length(parms))
    parms = rep(parms, reps)
    dfnames = paste0(fname,".", names(parms), ".", 1:length(parms), ".csv")
    dpnames = paste0(fname,".", names(parms), ".", 1:length(parms), ".parms.rds")
    a_ply(1:length(parms), 1, function(z) {
      saveRDS(parms[[z]], dpnames[[z]])
      run_sodi_single(parms=parms[[z]], init=init, progress=FALSE,
                      filename=dfnames[z])
    }, .progress = ifelse(progress, "time", "none"), .parallel=parallel)
  } else {
    saveRDS(parms, paste0(fname, ".parms.rds"))
    run_sodi_single(parms, init=init, progress=progress, filename=paste0(fname, ".csv"))
  }  
  if(load) {
    sodi <- load_sodi_files(dir = dir, name=fname)
  }
  setwd(wd)
  END = Sys.time()
  if(progress) cat("\nCompleted in:", END-START, "s\n")
  return(sodi)
}

#' @export
#' @import stringr
#' @importFrom Hmisc escapeRegex
load_sodi_files <- function(dir=dir, name=basename(dir)) {
  wd = getwd()
  setwd(dir)
  filenames = list.files(pattern=paste0(escapeRegex(name), ".*$"))
  dmatches <- str_match(filenames, "^([^\\.]*\\.?(\\w*)\\.?(\\d*))\\.csv$")
  dmatches = dmatches[complete.cases(dmatches),]
  pmatches <- str_match(filenames, "^(.*)\\.parms\\.rds$")
  pmatches = pmatches[complete.cases(pmatches),]
  if(!is.matrix(dmatches)) dmatches <- matrix(dmatches, nrow=1)
  if(!is.matrix(pmatches)) dmatches <- matrix(pmatches, nrow=1)
  sodi <- alply(dmatches[,2], 1, function(z) load_sodi_file_single(z, runs=(nrow(dmatches)!=1)))
  if(length(sodi)==1) {
    sodi = sodi[[1]]
    class(sodi) <- c(class(sodi), "sodi")
  } else {
    parmlist = llply(sodi, function(z) attr(z, "parms"))
    names(parmlist) = dmatches[,3]
    sodi <- rbindlist(sodi)
    attr(sodi, "parms") = parmlist[-which(duplicated(dmatches[,3]))]
    attr(sodi, "reps") = table(dmatches[,3])
    class(sodi) <- c(class(sodi), "sodi", "multisodi")

  }
  setwd(wd)
  return(sodi)
}


load_sodi_file_single <- function(name, runs=TRUE) {
  dfile = list.files(pattern=paste0("^", escapeRegex(name),  "\\.csv$"))
  pfile = list.files(pattern=paste0("^", escapeRegex(name),  "\\.parms\\.rds$"))
  if(length(dfile) > 1 || length(pfile) > 1) stop("Ambiguous filenames")
  colclasses=c("numeric", "integer", "numeric", "numeric", "integer", "integer")
  colnames=c("Time", "ID", "X", "Y", "Stage", "Infections")
  sodi <- fread(dfile, colClasses=colclasses)
  setnames(sodi, colnames)
  parms <- readRDS(pfile)
  sp = rep(1:length(parms$sp_stages), parms$sp_stages)
  ss = c(0, rep(which(diff(c(0, sp))==1), parms$sp_stages))
  sodi[, Species := parms$sp_names[sp[Stage[1]]], by=Stage]
  sodi[, Stage := Stage[1] - ss[Stage[1] + 1] + 1, by=Stage]
  if(runs) {
    match = str_match(name, "^.+\\.(\\w+)\\.(\\d+)$")
    sodi[, Parms := match[2]]
    sodi[, Run := as.numeric(match[3])]
  }
  attr(sodi, "parms") <- parms
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
  
  run_sodi_rcpp(init, parms_mod, progress, filename)
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
  