#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <memory>
#include <fstream>
#include <string>
#include <ctime>

#include "virus.hpp"
#include "host.hpp"
#include "hostpopulation.hpp"
#include "viruspopulation.hpp"

using namespace std;

int main(int argc, char *argv[]){
  int start = clock();
  srand(time(NULL));
  //for(int i = 0; i < 50; ++i){
  string filename;
  //filename = "output" + to_string(i) + ".csv";
  ofstream output ("output.csv");
  //  ofstream output(filename);
  ofstream voutput, voutput2;
  int day = 1;
  int final_day = 200;
  
  HostPopulation* hpop = new HostPopulation(90000,100,100000-90000-100,0,1.5,1.0/(40.0*365.0),1.0/25.0,0.333);
  while(day <= final_day){
    hpop->stepForward(day);
    hpop->printStatus();
    cout << endl;
    
    output << hpop->countSusceptibles() << "," << hpop->countInfecteds() << "," << hpop->countRecovereds() << endl;
    
    day++;
  }
  hpop->writeViruses(voutput, "voutput.csv");
  hpop->virusPairwiseMatrix(voutput2, "voutput2.csv",1000);
  delete hpop;
  output.close();

  int stop = clock();
  cout << "Time elapsed: " << (stop-start)/double(CLOCKS_PER_SEC) << " Seconds" << endl;
  //}
  return 0;
}
