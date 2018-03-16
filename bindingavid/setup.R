<<<<<<< HEAD
# preprocess                                                                                                                                                                                                                                                                                                                    # setup environment
rm(list = ls())
library("Rcpp")
library("reshape")

## install driftSim
#  from local disk
#setwd("F:/Documents/GitHub/antigenic_drift/CppInfectionSim/driftSim")
#devtools::load_all()
## from Github
setwd("F:/Documents/GitHub/antigenic_drift/bindingavid")
devtools::install_github("jameshay218/antigenic_drift/CppInfectionSim/driftSim")
library(driftSim)
packageVersion('driftSim')

## dependent files for batch analysis
myhome <- "F:/Documents/GitHub/antigenic_drift/bindingavid"
# Generate starting conditions
source(paste(c(myhome,"/cluster_submission.R"),collapse = ''))
source(paste(c(myhome,"/plot_SIR.R"),collapse = ''))
source(paste(c(myhome,"/saveoutput.R"),collapse = ''))

# start batch analysis with different scenarios
# (1)single_fixed_low, (2)single_fixed_high
# (3)single_adaptive_low, (4)single_adaptive_high, (5)multiple peaks
 source(paste(c(myhome,"/jobs/bindingavid_jobs1_pc.R"),collapse = '')) #(1)
# source(paste(c(myhome,"/jobs/bindingavid_jobs2_pc.R"),collapse = '')) #(2)
# source(paste(c(myhome,"/jobs/bindingavid_jobs3_pc.R"),collapse = '')) #(3)
# source(paste(c(myhome,"/jobs/bindingavid_jobs4_pc.R"),collapse = '')) #(4)
# source(paste(c(myhome,"/jobs/bindingavid_jobs6_2_pc.R"),collapse = '')) #(5)
=======
# preprocess                                                                                                                                                                                                                                                                                                                    # setup environment
rm(list = ls())
library("Rcpp")
library("reshape")

## install driftSim
#  from local disk
#setwd("F:/Documents/GitHub/antigenic_drift/CppInfectionSim/driftSim")
#devtools::load_all()
## from Github
setwd("F:/Documents/GitHub/antigenic_drift/bindingavid")
devtools::install_github("jameshay218/antigenic_drift/CppInfectionSim/driftSim")
library(driftSim)
packageVersion('driftSim')

## dependent files for batch analysis
myhome <- "F:/Documents/GitHub/antigenic_drift/bindingavid"
# Generate starting conditions
source(paste(c(myhome,"/cluster_submission.R"),collapse = ''))
source(paste(c(myhome,"/plot_SIR.R"),collapse = ''))
source(paste(c(myhome,"/saveoutput.R"),collapse = ''))

# start batch analysis with different scenarios
# (1)single_fixed_low, (2)single_fixed_high
# (3)single_adaptive_low, (4)single_adaptive_high
# (4.1) single_adaptive_high + different binding avidity 
# (6)multiple peaks
source(paste(c(myhome,"/jobs/bindingavid_jobs1_pc.R"),collapse = '')) #(1)
# source(paste(c(myhome,"/jobs/bindingavid_jobs2_pc.R"),collapse = '')) #(2)
# source(paste(c(myhome,"/jobs/bindingavid_jobs3_pc.R"),collapse = '')) #(3)
# source(paste(c(myhome,"/jobs/bindingavid_jobs4_pc.R"),collapse = '')) #(4)
# source(paste(c(myhome,"/jobs/bindingavid_jobs4_allbind_pc.R"),collapse = '')) #(4.1)
# source(paste(c(myhome,"/jobs/bindingavid_jobs6_2_pc.R"),collapse = '')) #(6)
>>>>>>> 20279a449d1ab5922b984792b494dae6bd7d6198
