## Read any files
#' runName: string for the name of the job (outputs files will be saved with this name)
#' runNo: number indicating which run this is
#' input_pars: a named vector of all the parameters needed for the simulation, as in input_params.csv"
#' deltaVMat: the full deltaV matrix as a matrix
#' flags: vector of booleans which are flags to certain settings in the run_simulation function
#' inputK: vector of initial Ks to create a starting distribution from
#' dur: duration of simulation in days
#' version: the version of the model to run as discussed previously
#' callback: this is to do with giving progress information back from C code, but just leave this as NULL
#' VERBOSE: bool indicating if we want console outputs
run_all <- function(runName,runNo,input_pars, deltaVMat,flags,inputK, dur,version, callback,VERBOSE){
####################
    ## Setup some pars
####################
  print(input_pars)
    
  print("Setup")
  
  if(!file.exists(paste("outputs/",runName,sep=""))) dir.create(paste("outputs/",runName,sep=""))
    filename1 <- paste("outputs/",runName,"/scenario_SIR_",runNo,".csv",sep="")
    filename2 <- paste("outputs/",runName,"/voutput1_",runNo,".csv",sep="")
    filename3 <- paste("outputs/",runName,"/voutput2_",runNo,".csv",sep="")
    filename4 <- paste("outputs/",runName,"/hosts_",runNo,".csv",sep="")
    filename5 <- paste("outputs/",runName,"/hostKs_",runNo,".csv",sep="")
    filename6 <- paste("outputs/",runName,"/SIR_",runNo,".png",sep="")
    filenames <- c(filename1, filename2, filename3, filename4, filename5)
    
####################

####################
    ## Read parameters
####################
    print("Read parameters")
    print(runName)
    all_pars <- input_pars[input_pars$runName == runName, c("names","values")]
    print(all_pars)
    print(runName)
    
    pars <- all_pars$values
    names(pars) <- all_pars$names

    hostpar_names <- c("s0","i0","r0","c","mu","w","g","iniBind","meanBoost","iniDist","saveFreq","maxTitre")
    viruspar_names <- c("p","r","q","a","b","n","v","probMut","expDist","kc","VtoD")

    hostpars <- pars[hostpar_names]
    viruspars <- pars[viruspar_names]
####################

####################
    ## Generate starting conditions
####################
    print("Generate starting conditions")
    #print(inputK)
    ## Read in the csv file and get the hostK distribution before calling the function
    hostFile <- inputFile
    hostDat <- read.csv(hostFile)
    print(hostFile)
    hostKs <- hostDat$hostK
    #print(hostKs)
    N <-hostpars["s0"]+hostpars["i0"]+hostpars["r0"]
    print(N)
    print(N)
    iniK <- generateHostKDist_2(hostFile,N)

#################### 
    ## Run simulation
####################
    print("run simulation")
    #readline()
    #run_simulation(flags=flags,hostpars=hostpars,viruspars=viruspars,deltaVMat=deltaVMat,iniKs=iniK,start=0,end=dur,input_k=unname(as.matrix(inputK)),output_files=filenames,VERBOSE=VERBOSE, scenario=version,callback=callback)
    save(flags, hostpars, viruspars, deltaVMat, iniK, dur, filenames, VERBOSE, version, callback, file = "runtime_objects.RData")
    run_simulation(flags=flags,hostpars=hostpars,viruspars=viruspars,deltaVMat=deltaVMat,iniKs=iniK,start=0,end=dur,output_files=filenames,VERBOSE=VERBOSE, scenario=version,callback=callback)

    
#####################
    ## SIR plot
#####################
    print(paste("plot ", filename1));
    dat <- read.csv(filename1)
    
    to.png <- function(expr, filename, ..., verbose=TRUE) {
      if ( verbose )
        cat(sprintf("Creating %s\n", filename))
      png(filename, ...)
      on.exit(dev.off())
      eval.parent(substitute(expr))
    }
    print("finish ploting");
    to.png(plot(plot_SIR(filename1,N)),filename6)
    ####################
    
    
    
    
  
    return(TRUE)
}

## Add this to the source file so that you can access it ie. cluster_submission.R
generateHostKDist_2<- function(hostKs, N){
  countHostK <- count(hostKs)
  freqs <- countHostK$freq/sum(countHostK$freq)
  cumSumK <- cumsum(freqs)
  startingKs <- generateKSamples(cumSumK, N)
  startingKs <- startingKs + 1
return(countHostK$x[startingKs])}

