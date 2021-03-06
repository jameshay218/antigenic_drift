library(Rcpp)
library(ggplot2)

setup_env <- function(){
    Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
    sourceCpp("~/Documents/antigenic_drift/antigenic_drift/CppInfectionSim/src/rcpp_interface.cpp", rebuild=TRUE)
}

run_simulation <- function(
                           flags=c(1,0,0,0),
                           hostpars=c(90000,100,100000-90000-100,1.5,1/(40*365),1/25,0.333,0.8),
                           viruspars=c(3, 50, 1, 0.7, 3, 2, 2, 0.1, 1, 0.5, 1000),
                           start=0,
                           end=365,
                           output_files = c("SIR_output.csv","voutput_1.csv","voutput_2.csv"),
                           VERBOSE=TRUE,
                           scenario){
    result = tryCatch({
        run_simulation_cpp(flags,hostpars,viruspars,start,end,output_files,VERBOSE,scenario=1)
    }, warning=function(w){
        print("Warning: ")
        print(w)
    }, error = function(e){
        print("Error - have you run setup_env?")
        print(e)
    }, finally = {
    }
    )
    
}
