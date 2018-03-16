library(driftSim)

setwd("~/tmp/scripts/")

runs <- 100

SIR_flag <- 1
voutput1_flag <- 0 #' Flag to save virus information for Sean's phylogenetic tree
voutput2_flag <- 0 #' Flag to save pairwise distance matrix
time_flag <- 0 #' Flag to record time taken for simulation
VERBOSE <- 0 #' Outputs in simulation
save_state <- 0 #' Flag to save the final state of the simulation
input_flagA <- 0 #' Flag to use specified file as input for simulation
input_flagB <- 1 #' Flag to use specified file as input for simulation
save_k <- 0

flags <- c(SIR_flag, voutput1_flag, voutput2_flag, time_flag, save_state, input_flagA, input_flagB, save_k)

S0 <- 999900
I0 <- 100
R0 <- 0

contactRate <- 0.7
mu <- 1/(70*365)
wane <- 1/(25)
gamma <- 1/3.3
##################
iniBindA <- 0.45
iniBindB <- 0.6
##################
meanBoost <- 6
iniDist <- 2
saveFreq <- 5

duration <- 1000


deltaVMat <- unname(as.matrix(read.csv("deltaVMat.csv",header=FALSE)))

p <- 4
q <- 1
r <- 70
b <- 3
a <- 0.7
n <- 4
v <- 1

probMut <- 0
expDist <- 1
###################
kc <- 0.3
###################
VtoD <- 0.1

viruspars <- c(p,r,q,a,b,n,v,probMut,expDist,kc,VtoD)

inputFiles <- c("iniK.csv")
iniK <- read.csv(inputFiles[1],header=FALSE)[,1]

for(i in 1:runs){
    hostpars <- c(S0,I0, R0,contactRate,mu,wane,gamma,iniBindA, meanBoost, iniDist,saveFreq)
    print(paste("Run number: ", i,sep=""))
    filename1 <- paste("out/SIRfixed_low_",i,".csv",sep="")
    filename2 <- paste("voutput1_",i,".csv",sep="")
    filename3 <- paste("voutput2_",i,".csv",sep="")
    filename4 <- paste("hosts_",i,".csv",sep="")
    filename5 <- paste("hostKs_",i,".csv",sep="")
    filenames <- c(filename1, filename2, filename3, filename4, filename5)
    print("Fixed V low...")
    y <- run_simulation(flags,hostpars,viruspars,deltaVMat,iniK,0,duration,inputFiles,filenames,VERBOSE, 1,NULL)
    filenames[1] <- paste("out/SIRadaptive_low_",i,".csv",sep="")
    print("Adaptive V low...")
    y <- run_simulation(flags,hostpars,viruspars,deltaVMat,iniK,0,duration,inputFiles,filenames,VERBOSE, 2,NULL)
    
    hostpars <- c(S0,I0, R0,contactRate,mu,wane,gamma,iniBindB, meanBoost, iniDist,saveFreq)
     print("Fixed V high...")
    filenames[1] <- paste("out/SIRfixed_high_",i,".csv",sep="")
    y <- run_simulation(flags,hostpars,viruspars,deltaVMat,iniK,0,duration,inputFiles,filenames,VERBOSE, 1,NULL)
    filenames[1] <- paste("out/SIRadaptive_high_",i,".csv",sep="")
    print("Adaptive V high...")
    y <- run_simulation(flags,hostpars,viruspars,deltaVMat,iniK,0,duration,inputFiles,filenames,VERBOSE, 2,NULL)
    
    
}
