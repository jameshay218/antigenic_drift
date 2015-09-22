#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <memory>

#include "virus.hpp"
#include "host.hpp"
#include "hostpopulation.hpp"
#include "viruspopulation.hpp"

using namespace std;

int main(int argc, char *argv[]){
  int day = 1;
  int final_day = 1000;
  HostPopulation* hpop = new HostPopulation(1000000,10,0,0,0.05,0.001,0.01,0.01);

  while(day <= final_day){
    hpop->stepForward(day);
    hpop->printStatus();
      
    day++;
  }
  return 0;
}
