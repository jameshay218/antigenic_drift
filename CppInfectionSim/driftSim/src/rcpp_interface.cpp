//rcpp_functions.cpp
//[[Rcpp::plugins(cpp11)]]

#include <iostream>
#include <vector>
#include <time.h>
#include <memory>
#include <fstream>
#include <string>
#include <ctime>
#include <Rmath.h>
#include <Rcpp.h>
#include <vector>

#include "virus.hpp"
#include "host.hpp"
#include "hostpopulation.hpp"

using namespace std;

//[[Rcpp::export]]
int run_simulation_cpp(Rcpp::IntegerVector flags, Rcpp::NumericVector hostPopn, Rcpp::NumericVector virusPars, int day, int final_day,std::vector<std::string> output_files, bool VERBOSE,int scenario, SEXP callback){
  bool save_SIR = flags[0];
  bool save_viruses = flags[1];
  bool save_pairwise_viruses = flags[2];
  bool use_time = flags[3];

  double S0 = hostPopn[0];
  double I0 = hostPopn[1];
  double R0 = hostPopn[2];
  double start_day = day;
  double contactRate = hostPopn[3];
  double mu = hostPopn[4];
  double wane = hostPopn[5];
  double gamma = hostPopn[6];
  double iniBinding = hostPopn[7];

  double p = virusPars[0];
  double r = virusPars[1];
  double q = virusPars[2];
  double a = virusPars[3];
  double b = virusPars[4];
  double n = virusPars[5];
  double v = virusPars[6];
  double probMut = virusPars[7];
  double expDist = virusPars[8];
  double kc = virusPars[9];
  double VtoD = virusPars[10];

  int start = day;
  
  const bool has_callback = callback != R_NilValue;
  if (has_callback) {
    callback = Rcpp::as<Rcpp::Function>(callback);
  }

  Virus::set_p(p);
  Virus::set_r(r);
  Virus::set_q(q);
  Virus::set_a(a);
  Virus::set_b(b);
  Virus::set_n(n);
  Virus::set_v(v);
  Virus::set_prob_mut(probMut);
  Virus::set_exp_dist(expDist);
  Virus::set_kc(kc);
  Virus::set_VtoD(VtoD);
  Virus::set_scenario(scenario);


  if(use_time) int start = clock();
  string filename;
 ofstream output (output_files[0]);
 ofstream voutput;
 ofstream voutput2;
 if(VERBOSE){
   Rcpp::Rcout << "Params: " << endl;
   Rcpp::Rcout << "S0: " << S0 << endl;
   Rcpp::Rcout << "I0: " << I0 << endl;
   Rcpp::Rcout << "R0: " << I0 << endl;
   Rcpp::Rcout << "contact: " << contactRate << endl;
   Rcpp::Rcout << "mu: " << mu << endl;
   Rcpp::Rcout << "wane: " << wane << endl;
   Rcpp::Rcout << "gamma: " << gamma << endl;
   Rcpp::Rcout << "iniBinding: " << iniBinding << endl;
   Rcpp::Rcout << "Scenario: " << scenario << endl;
 }

 HostPopulation* hpop = new HostPopulation(
					   S0,
					   I0,
					   R0,
					   start_day,
					   contactRate,
					   mu,
					   wane,
					   gamma,
					   iniBinding
					   );
 while(day <= final_day){
   hpop->stepForward(day);
   if(VERBOSE){
     hpop->printStatus();
     Rcpp::Rcout << endl;
   }
    
   if(save_SIR) output << hpop->countSusceptibles() << "," << hpop->countInfecteds() << "," << hpop->countRecovereds() << endl;
   if (has_callback) {
     Rcpp::IntegerVector tmp = Rcpp::IntegerVector::create(
							   Rcpp::_["day"]=hpop->getDay(),
							   Rcpp::_["susceptibles"]=hpop->countSusceptibles(),
							   Rcpp::_["infecteds"]=hpop->countInfecteds(),
							   Rcpp::_["recovereds"]=hpop->countRecovereds());
     Rcpp::as<Rcpp::Function>(callback)(tmp);
   }
   day++;
 }

 if(save_viruses) hpop->writeViruses(voutput, output_files[1]);
 if(save_pairwise_viruses) hpop->virusPairwiseMatrix(voutput2, output_files[2],1000);
 delete hpop;
 if(save_SIR) output.close();
  
 if(use_time){
   int stop = clock();
   Rcpp::Rcout << "Time elapsed: " << (stop-start)/double(CLOCKS_PER_SEC) << " Seconds" << endl;
 }
 return 0;
}
