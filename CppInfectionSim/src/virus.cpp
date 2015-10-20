#include "virus.hpp"
#include "host.hpp"
#include<iostream>
#include<math.h>
#include<vector>

using namespace std;

int Virus::v_IDgenerator;
double Virus::_p;
double Virus::_r;
double Virus::_q;
double Virus::_a;
double Virus::_b;
double Virus::_n;
double Virus::_v;
double Virus::_k;


Virus::Virus(Virus* _parent, double _bindavid, double _distance, Host* _host, int _t){
  id = v_IDgenerator++;
  birth = _t;
  infectionK = 0;
  bindingavid = _bindavid;
  parent = _parent;
  level = parent->level + 1;
  distanceToParent = _distance;
  host = _host;
}

Virus::Virus(int _level, Virus* _parent, double _bindavid, double _distance, Host* _host, int _t){
  id = v_IDgenerator++;
  birth = _t;
  infectionK = 0;
  parent = _parent;
  bindingavid = _bindavid;
  level = _level;
  distanceToParent = _distance;
  host = _host;
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

int Virus::getBirth(){
  return birth;
}

int Virus::getDeath(){
  return death;
}

int Virus::getK(){
  return infectionK;
}

void Virus::updateK(int _k){
  infectionK = _k;
}

void Virus::kill(int cur_t){
  death = cur_t;
}

double Virus::calculateRho(){
  return(1 - pow((1/(probSurvival()*probReplication())),(-_n*_v)));
}

void Virus::mutate(){
}

double Virus::probSurvival(){
  double sum = 0;
  int size  = host->getInfectionHistory().size();
  for(int i = 0; i < size; ++i){
    sum += _r*_k - getAntigenicDistance(this, host->getInfectionHistory()[i]);
  }
  return(pow((1 - exp(-_p*(bindingavid + _q))), sum)); 
}

double Virus::probReplication(){
  return(exp(-_a*(pow(bindingavid,_b))));
}
	 

double Virus::getAntigenicDistance(Virus* A, Virus* B){
    int levelA = A->getLevel();
    int levelB = B->getLevel();
    int tmp;
    Virus* tmpA;
    Virus* tmpB;
    double totalDistance = 0;

    // Swap so level A is lower (ie. level is a higher int) than level B
    if(levelB > levelA){
      tmp = levelA;
      levelA = levelB;
      levelB = tmp;
      tmpA = A;
      A = B;
      B = tmpA;
    }

    // Starting at the lower level, backtrack until at same level as other virus
    tmpA = A;
    while(tmpA->getLevel() > levelB){
      totalDistance += tmpA->getDistance();
      tmpA = tmpA->getParent();
    }

    // If tmpA is now at the other virus, then B is an ancestor of A and we can return
    if(tmpA == B){ return totalDistance; }

    // Otherwise, backtrack from both at the same time until we find a cmomon ancestor. If either reaches NULL, then they are not related.
    tmpB = B;
    totalDistance += tmpA->getDistance() + tmpB->getDistance();

    // Need to add a check for NULL parent (ie. different trees, so infinite distance)
    while((tmpA->getParent() != NULL && tmpB->getParent() != NULL) && tmpA->getParent() != tmpB->getParent()){
      tmpA = tmpA->getParent();
      tmpB = tmpB->getParent();
      totalDistance += tmpA->getDistance() + tmpB->getDistance();
    }


    return totalDistance;
  };
