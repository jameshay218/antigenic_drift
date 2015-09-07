#include "virus.hpp"
#include<iostream>
using namespace std;

Virus::Virus(int _id, Virus* _parent, double _bindavid, double _distance, Host* _host, VirusPopulation* popn){
  id = _id;
  bindingavid = _bindavid;
  parent = _parent;
  level = parent->level + 1;
  distanceToParent = _distance;
}


Virus::Virus(int _id, int _level, Virus* _parent, double _bindavid, double _distance, Host* _host, VirusPopulation* popn){
  id = _id;
  parent = _parent;
  bindingavid = _bindavid;
  level = _level;
  distanceToParent = _distance;
}


int Virus::getId(){
  return id;
}

Virus* Virus::getParent(){
  return parent;
}

int Virus::getLevel(){
  return level;
}

double Virus::getDistance(){
  return distanceToParent;
}

double Virus::getBindingAvid(){
  return bindingavid;
}

double Virus::calculateRho(){
  return 0;
}

void Virus::mutate(){
}

double Virus::probSurvival(){
  return 0;
}

double Virus::probReplication(){
  return 0;
}
