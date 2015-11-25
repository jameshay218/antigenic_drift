#include "virus.hpp"
#include "host.hpp"
#include<iostream>
#include<math.h>
#include<vector>

using namespace std;

int Virus::v_IDgenerator;
double Virus::_p=1;
double Virus::_r=1;
double Virus::_q=1;
double Virus::_a=0.7;
double Virus::_b=3;
double Virus::_n=2;
double Virus::_v=1;
double Virus::_k=1;


Virus::Virus(Virus* _parent, Host* _host, int _t){
  id = v_IDgenerator++;
  birth = _t;
  infectionK = 0;
  bindingavid = bindingavid_ini = _parent->getBindingAvid();
  parent = _parent;
  level = parent->level + 1;
  distanceToParent = 0;
  host = _host;
}

Virus::Virus(int _level, Virus* _parent, double _bindingavid, double _distance, Host* _host, int _t){
  id = v_IDgenerator++;
  birth = _t;
  infectionK = 0;
  parent = _parent;
  bindingavid = bindingavid_ini = _bindingavid;
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

double Virus::getIniBindingAvid(){
  return bindingavid_ini;
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

double Virus::calculateRho(Host* _host){
  return(1 - pow((1/(_n*probSurvival(_host)*probReplication())),_v));
}

void Virus::mutate(){
  double prob_mut = 0.00;
  double tmp = ((double) rand() / (RAND_MAX));
  std::default_random_engine exp_gen;
  std::exponential_distribution<double> dist(1);
  double change = dist(exp_gen);

  if(tmp <= prob_mut){
    distanceToParent += change;
  }
}

double Virus::probSurvival(Host* _host){
  double sum = 0;
  int size  = _host->getInfectionHistory().size();
  for(int i = 0; i < size; ++i){
    sum += _r*_k - getAntigenicDistance(this, _host->getInfectionHistory()[i]);
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

    // Otherwise, backtrack from both at the same time until we find a common ancestor. If either reaches NULL, then they are not related.
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
