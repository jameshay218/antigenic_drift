# setup environment
rm(list = ls())
#setwd("F:/Documents/GitHub/antigenic_drift/CppInfectionSim/driftSim")
#devtools::load_all()
library("reshape")
# go to antigenic_drift/bindingavid subfolder
setwd("F:/Documents/GitHub/antigenic_drift/bindingavid")
devtools::install_github("jameshay218/antigenic_drift/CppInfectionSim/driftSim")
library(driftSim)


# start batch analysis
# update the source folder
myhome <- "F:/Documents/GitHub/antigenic_drift/bindingavid"
source(paste(c(myhome,"/cluster_submission.R"),collapse = ''))
source(paste(c(myhome,"/plot_SIR.R"),collapse = ''))
source(paste(c(myhome,"/bindingavid_jobs.R"),collapse = ''))

