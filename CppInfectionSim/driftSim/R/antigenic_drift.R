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
                           flags=c(1,0,0,0),
                           hostpars=c(90000,100,100000-90000-100,1.5,1/(40*365),1/25,0.333,0.8),
                           viruspars=c(3, 50, 1, 0.7, 3, 2, 2, 0.1, 1, 0.5, 1000),
                           start=0,
                           end=365,
                           output_files = c("SIR_output.csv","voutput_1.csv","voutput_2.csv"),
                           VERBOSE=TRUE,
                           scenario=1,
                           callback=NULL){
    print("Here")
    run_simulation_cpp(flags,hostpars,viruspars,start,end,output_files,VERBOSE,scenario, callback)
}
