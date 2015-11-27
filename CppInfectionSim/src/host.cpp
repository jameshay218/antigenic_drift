#include "host.hpp"

using namespace std;

Host::Host(){
  state = Susceptible;
  currentInfection = NULL;
  hostK = 0;
}

Host::Host(State _state, HostPopulation* _popn){
  state = _state;
  currentInfection = NULL;
  popn = _popn;
  hostK = 0;
}

Host::~Host(){
  int j = infectionHistory.size();
  for(int i = 0; i < j;++i){
    delete infectionHistory[i];
  }
}

void Host::infect(Virus* newInfection, int cur_t){
  state = Infected;

  newInfection->updateK(infectionHistory.size()+1);
  if(currentInfection != NULL){
    currentInfection->kill(cur_t);
    infectionHistory.push_back(currentInfection);
  }
  currentInfection = newInfection;
}

void Host::recover(int cur_t){
  state = Recovered;

  poisson_distribution<int> poisson(5);
  int boost = poisson(popn->generator);
  hostK += boost;

  if(currentInfection != NULL){
    currentInfection->kill(cur_t);
    infectionHistory.push_back(currentInfection);
  }
  currentInfection = NULL;
}

Virus* Host::getCurrentVirus(){
  return(currentInfection);
}

std::vector<Virus*> Host::getInfectionHistory(){
  return(infectionHistory);
}

State Host::getState(){
  return(state);
}

void Host::die(int cur_t){
  state = Dead;
  if(currentInfection != NULL){
    currentInfection->kill(cur_t);
    infectionHistory.push_back(currentInfection);
  }
  currentInfection = NULL;
}

void Host::wane(){
  state = Susceptible;
}

double Host::calculateBeta(){
  return(popn->getContactRate()*currentInfection->calculateRho(this));
}

double Host::get_hostK(){
  return(hostK);
}

default_random_engine Host::get_generator(){
  return(popn->generator);
}
