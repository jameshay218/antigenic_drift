#include "host.hpp"

using namespace std;

Host::Host(){
  state = Susceptible;
  currentInfection = NULL;
}

Host::Host(State _state, HostPopulation* _popn){
  state = _state;
  currentInfection = NULL;
  popn = _popn;
}

void Host::infect(Virus* newInfection){
  state = Infected;
  if(currentInfection != NULL){
    infectionHistory.push_back(currentInfection);
  }
  currentInfection = newInfection;
}

void Host::recover(){
  state = Recovered;
  if(currentInfection != NULL){
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

void Host::die(){
  state = Dead;
  if(currentInfection != NULL){
    infectionHistory.push_back(currentInfection);
  }
  currentInfection = NULL;
}

void Host::wane(){
  state = Susceptible;
}

double Host::calculateBeta(){
  return(popn->getContactRate()*currentInfection->calculateRho());
}
