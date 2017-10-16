#update myhome according to your system
#myhome <- "~/Documents/Binding Avidity/antigenic_drift/bindingavid/"

#source("~/net/home/bindingavid/scripts/cluster_setup.R")
#source("scripts/cluster_setup.R")
#source("~/net/home/bindingavid/scripts/cluster_submission.R")
setwd(myhome)
source(paste(c(myhome,"/cluster_submission.R"),collapse = ''))
library("plyr")
library(ggplot2)
library(reshape2)
source("plot_SIR.R")
#devtools::load_all("~/Documents/Binding Avidity/antigenic_drift/CppInfectionSim/driftSim/")

#setwd("~/net/home/bindingavid")
#setwd("F:/Documents/GitHub/antigenic_drift/bindingavid")
#myhome <- "F:/Documents/GitHub/antigenic_drift/bindingavid"


param.types = c("single_fixed_low","single_fixed_high","single_adaptive_low","single_adaptive_high","multiple_fixed","multiple_adaptive_test2")      
param.filenames = c("input_params_se01_SFL.csv","input_params_se02_SFH.csv","input_params_se03_SAL.csv","input_params_se04_SAH.csv","input_params_se05_MF.csv","input_params_se06_2_MA.csv")
params.inputFile = c("hosts_1_ini.csv","hosts_1_ini.csv","hosts_1_ini.csv","hosts_1_ini.csv","hosts_3_ini2.csv","hosts_1_ini.csv")

#For single epidemic
#runName <- "single_fixed_low"
#runName <- "single_fixed_high"
#runName <- "single_adaptive_low"
runName <- "single_adaptive_high"
#runName <- "multiple_fixed"
#runName <- "multiple_adaptive_test2"

flag_run <- match(runName,param.types)

version <- 0 # it represents the scenario case
# Scenario 1: Random drift; fixed binding avidity
# Scenario 2: No drift; adaptive binding avidity; adaptive antigenic change
# Scenario 3: Random drift; adaptive binding avidity; adaptive antigenic change
# Scenario 4: Random drift; adaptive binding avidity; no adaptive antigenic change
# Only use sceniario 1 and 3


if (flag_run == 1){ #"single_fixed_low" 
  version <- 1
  inputFile <- paste(c(myhome, "/inputs/", params.inputFile[flag_run]),collapse = '')
  input_pars <- read.csv(paste(c(myhome,"/inputs/", param.filenames[flag_run]),collapse = ''),stringsAsFactors=FALSE)
}

if (flag_run == 2){ #"single_fixed_high" 
  version <- 1
  print("single_fixed_high")
  inputFile <- paste(c(myhome, "/inputs/", params.inputFile[flag_run]),collapse = '')
  input_pars <- read.csv(paste(c(myhome,"/inputs/", param.filenames[flag_run]),collapse = ''),stringsAsFactors=FALSE)
}

if (flag_run == 3 | flag_run == 4 | flag_run == 5 | flag_run == 6 ){ 
  version <- 3
  inputFile <- paste(c(myhome, "/inputs/", params.inputFile[flag_run]),collapse = '')
  input_pars <- read.csv(paste(c(myhome,"/inputs/", param.filenames[flag_run]),collapse = ''),stringsAsFactors=FALSE)
}

#if (flag_run == 4){ #"single_adaptive_high" 
#  version <- 3
#  inputFile <- paste(c(myhome, "/inputs/", params.inputFile[flag_run]),collapse = '')
#  input_pars <- read.csv(paste(c(myhome,"/inputs/", param.filenames[flag_run]),collapse = ''),stringsAsFactors=FALSE)
#}



deltaVMat <- unname(as.matrix(read.csv(paste(c(myhome,"/inputs/deltaVMat.csv"),collapse = ''),sep=",",stringsAsFactors=FALSE),header=FALSE))
inputK <- read.csv(inputFile,stringsAsFactors=FALSE)



dur <- 2000
callback <- NULL

runs <- 200

#################
## Setup output flags
#################
VERBOSE <- FALSE
SIR_flag <- TRUE #' Flag to save SIR dynamics
voutput1_flag <- FALSE #' Flag to save virus information for Sean's phylogenetic tree
voutput2_flag <- FALSE #' Flag to save pairwise distance matrix
time_flag <- FALSE #' Flag to record time taken for simulation
VERBOSE <- FALSE #' Outputs in simulation
save_state <- FALSE #' Flag to save the final state of the simulation
input_flag_generated <- TRUE #' Flag to use specified file as input for simulation
input_flag_saved <- FALSE #' Flag to use specified file as input for simulation
save_k <- FALSE
flags <- c(SIR_flag, voutput1_flag, voutput2_flag, time_flag, save_state, input_flag_generated, input_flag_saved, save_k)
flags <- as.numeric(flags)
####################

combos <- expand.grid(runName=runName,runNo=1:runs,stringsAsFactors=FALSE)

if(!file.exists(paste(myhome, "/outputs/",runName,sep=""))) dir.create(paste(myhome, "/outputs/",runName,sep=""))

log <- ""
Kc = seq(0, 0.6, by=0.05)
#Kc <- c(0, 0.6, 0.05)
run_id <- 1
for(i in 1:length(Kc)) {
  print(Kc[i])
  input_pars[input_pars[,2] == "kc",3]<- Kc[i]
  #vtod_par <- input_pars[input_pars[,2] == "VtoD",3]
  #vtod <- (vtod_par*(0.05/Kc[i]))
  #input_pars[input_pars[,2] == "VtoD",3]<- vtod
  vtod <- input_pars[input_pars[,2] == "VtoD",3]
  expDist <- input_pars[input_pars[,2] == "expDist",3]
  cont <- input_pars[input_pars[,2] == "c",3]
  for (runNo in 1:runs) {
    
    run_all(runName=runName,runNo=run_id,input_pars=input_pars,deltaVMat=deltaVMat,
            flags=flags,inputK=inputK,dur=dur,version=version,callback=callback,VERBOSE=VERBOSE)
    log <- paste(log, "run", run_id, ", kc=", Kc[i], ", VtoD=", vtod, ", expDist=", expDist, ",cont=", cont, " \n", sep="")
    run_id <- run_id + 1
  }
  
}


#the following code is for cluster at IC
#submission <- queuer::enqueue_bulk(obj1,combos, "run_all",input_pars=input_pars,deltaVMat=deltaVMat,flags=flags,inputK=inputK,dur=dur,version=version,callback=callback,VERBOSE=VERBOSE,do.call=TRUE,timeout=0)


#####################
##create output sub folders
#####################
outputdir <- paste(myhome,"/outputs/",runName,"/",sep="")
write(log,paste(outputdir,file="log.txt",sep=""))
print('print log file.')
saveoutput(outputdir, runName, input_pars)



