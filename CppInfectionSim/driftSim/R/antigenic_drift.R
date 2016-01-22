##' moar
##' @title xx
##' @param flags xx
##' @param hostpars xx
##' @param viruspars xx
##' @param start xx
##' @param end xx
##' @param output_files xx
##' @param VERBOSE xx
##' @param scenario xx
##' @param callback xx
##' @export
##' @import Rcpp
##' @useDynLib driftSim
run_simulation <- function(
                           flags=c(1,0,0,0,0,0),
                           hostpars=c(90000,100,100000-90000-100,1.5,1/(40*365),1/25,0.333,0.8,10),
                           viruspars=c(3, 70, 1, 0.7, 3, 2, 2, 0.1, 1, 0.5, 1000),
                           deltaVMat = NULL,
                           start=0,
                           end=365,
                           input_files=c("hosts.csv","viruses.csv"),
                           output_files = c("SIR_output.csv","voutput_1.csv","voutput_2.csv","hosts.csv","viruses.csv"),
                           VERBOSE=TRUE,
                           scenario=1,
                           callback=NULL){
    run_simulation_cpp(flags,hostpars,viruspars, deltaVMat, start,end,input_files,output_files,VERBOSE,scenario, callback)
}
