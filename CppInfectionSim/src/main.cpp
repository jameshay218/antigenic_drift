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
  ofstream output ("output.csv");
  ofstream voutput;
  int day = 1;
  int final_day = 500;
  HostPopulation* hpop = new HostPopulation(100000,10,0,0,0.05,0.001,0.01,0.01);

  while(day <= final_day){
    hpop->stepForward(day);
    hpop->printStatus();
    
    output << hpop->countSusceptibles() << "," << hpop->countInfecteds() << "," << hpop->countRecovereds() << endl;
    
    day++;
  }
  hpop->writeViruses(voutput, "voutput.csv");
  output.close();

  return 0;
}
