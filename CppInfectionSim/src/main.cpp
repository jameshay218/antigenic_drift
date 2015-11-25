#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <memory>
#include <fstream>
#include <string>

#include "virus.hpp"
#include "host.hpp"
#include "hostpopulation.hpp"
#include "viruspopulation.hpp"

using namespace std;

int main(int argc, char *argv[]){
  srand(time(NULL));
  ofstream output ("output.csv");
  ofstream voutput;
  int day = 1;
  int final_day = 1000;
  HostPopulation* hpop = new HostPopulation(90000,70,9930,0,0.9,(1.0/(40.0*365.0)),0.04,0.263);

  while(day <= final_day){
    hpop->stepForward(day);
    hpop->printStatus();
    cout << endl;
    
    output << hpop->countSusceptibles() << "," << hpop->countInfecteds() << "," << hpop->countRecovereds() << endl;
    
    day++;
  }
  hpop->writeViruses(voutput, "voutput.csv");
  output.close();

  return 0;
}
