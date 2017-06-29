myhome <- "~/Documents/Binding Avidity/antigenic_drift/bindingavid/"
#source("~/net/home/bindingavid/scripts/cluster_setup.R")
#source("scripts/cluster_setup.R")
#source("~/net/home/bindingavid/scripts/cluster_submission.R")
setwd(myhome)
source(paste(c(myhome,"/cluster_submission.R"),collapse = ''))
library("plyr")
library(ggplot2)
library(reshape2)
source("plot_SIR.R")
devtools::load_all("~/Documents/Binding Avidity/antigenic_drift/CppInfectionSim/driftSim/")
#setwd("~/net/home/bindingavid")
#setwd("F:/Documents/GitHub/antigenic_drift/bindingavid")
#myhome <- "F:/Documents/GitHub/antigenic_drift/bindingavid"


param.types = c("single_fixed_low","single_fixed_high","single_adaptive","multiple_fixed","multiple_adaptive")      
param.filenames = c("input_params_se01_SFL.csv","input_params_se01_SFH.csv","input_params_SA.csv","input_params_MF.csv","input_params_MA.csv")
params.inputFile = c("hosts_1_ini.csv","hosts_1_ini.csv","hosts_1_ini.csv","hosts_3_ini2.csv","hosts_3_ini2.csv")

#For single epidemic
runName <- "single_fixed_low"
flag_run <- match(runName,param.types)
#if (flag_run == 1){
  inputFile <- paste(c(myhome, "inputs/", params.inputFile[flag_run]),collapse = '')
  input_pars <- read.csv(paste(c(myhome,"/inputs/", param.filenames[flag_run]),collapse = ''),stringsAsFactors=FALSE)
#}

#if (flag_run == 2){
#  inputFile <- paste(c(myhome, "/inputs/hosts_1_ini.csv"),collapse = '')
#  input_pars <- read.csv(paste(c(myhome,"/inputs/", param.filenames[flag_run]),collapse = ''),stringsAsFactors=FALSE)
#}

deltaVMat <- unname(as.matrix(read.csv(paste(c(myhome,"/inputs/deltaVMat.csv"),collapse = ''),sep=",",stringsAsFactors=FALSE),header=FALSE))
inputK <- read.csv(inputFile,stringsAsFactors=FALSE)

version <- 1
dur <- 800
callback <- NULL

runs <- 1

#################
## Setup output flags
#################
VERBOSE <- FALSE
SIR_flag <- TRUE #' Flag to save SIR dynamics
voutput1_flag <- TRUE #' Flag to save virus information for Sean's phylogenetic tree
voutput2_flag <- TRUE #' Flag to save pairwise distance matrix
time_flag <- FALSE #' Flag to record time taken for simulation
VERBOSE <- FALSE #' Outputs in simulation
save_state <- TRUE #' Flag to save the final state of the simulation
input_flag_generated <- TRUE #' Flag to use specified file as input for simulation
input_flag_saved <- FALSE #' Flag to use specified file as input for simulation
save_k <- TRUE
flags <- c(SIR_flag, voutput1_flag, voutput2_flag, time_flag, save_state, input_flag_generated, input_flag_saved, save_k)
flags <- as.numeric(flags)
####################

combos <- expand.grid(runName=runName,runNo=1:runs,stringsAsFactors=FALSE)

if(!file.exists(paste(myhome, "outputs/",runName,sep=""))) dir.create(paste(myhome, "outputs/",runName,sep=""))

for (runNo in 1:runs) {
  run_all(runName=runName,runNo=runNo,input_pars=input_pars,deltaVMat=deltaVMat,
          flags=flags,inputK=inputK,dur=dur,version=version,callback=callback,VERBOSE=VERBOSE)
}

#the following code is for cluster at IC
#submission <- queuer::enqueue_bulk(obj1,combos, "run_all",input_pars=input_pars,deltaVMat=deltaVMat,flags=flags,inputK=inputK,dur=dur,version=version,callback=callback,VERBOSE=VERBOSE,do.call=TRUE,timeout=0)


#####################
##create sub folders
#####################
print("create subfolders")
outdir <- paste("outputs/",runName,sep="")
subdir = list.dirs(path = outdir, full.names = TRUE, recursive = TRUE)
subfolders <- sapply(strsplit(subdir, "/"), "[", 3)
subfolders[is.na(subfolders)]<-0
v <- as.numeric(subfolders)
maxv <- max(v*1000)  ## to do 20170620
nextv = maxv + 1
newfolder <- ""
if (nextv >=10) {
  if (nextv >=100) {
    newfolder = paste(".", toString(nextv),sep = "")
  } else {
    newfolder = paste(".0", toString(nextv),sep = "")  
  }
} else {
  newfolder = paste(".00", toString(nextv),sep = "")  
}
newpath <- paste(outdir,"/",toString(newfolder),sep="")
dir.create(newpath)
print(newpath)
if (file.exists(newpath)) {
  filelist <- dir(outdir)
  for (fileID in 1:length(filelist)) {
    #print(filelist[fileID])
    file.rename(from=paste(outdir,"/",filelist[fileID],sep=""),to=paste(newpath,"/",filelist[fileID],sep=""))
  }
}



