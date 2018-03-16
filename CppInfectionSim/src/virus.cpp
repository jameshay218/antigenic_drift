#include "virus.hpp"
#include "host.hpp"
#include<iostream>
#include<math.h>
#include<vector>

using namespace std;

int Virus::v_IDgenerator = 1;
double Virus::_p=3.0;
double Virus::_r=50.0;
double Virus::_q=1.0;
double Virus::_a=0.7;
double Virus::_b=3.0;
double Virus::_n=2.0;
double Virus::_v=2.0;
double Virus::_prob_mut = 0.1;
double Virus::_exp_dist = 1;
double Virus::_kc = 0.5;
double Virus::_V_to_d = 1000;
int Virus::_scenario = 1;

void Virus::set_default(){
  _p=3.0;
  _r=50.0;
  _q=1.0;
  _a=0.7;
  _b=3.0;
  _n=2.0;
  _v=2.0;
  _prob_mut = 0.1;
  _exp_dist = 1;
  _kc = 0.5;
  _V_to_d = 1000;
  _scenario = 1;
}
void Virus::set_scenario(int _scen){
  _scenario = _scen;
}

void Virus::set_p(double new_p){
  _p = new_p;
}
void Virus::set_r(double new_r){
  _r = new_r;
}
void Virus::set_q(double new_q){
  _q = new_q;
}
void Virus::set_a(double new_a){
  _a = new_a;
}
void Virus::set_b(double new_b){
  _b = new_b;
}
void Virus::set_n(double new_n){
  _n = new_n;
}
void Virus::set_v(double new_v){
  _v = new_v;
}
void Virus::set_prob_mut(double new_probMut){
  _prob_mut = new_probMut;
}
void Virus::set_exp_dist(double new_exp){
  _exp_dist = new_exp;
}
void Virus::set_kc(double new_kc){
  _kc = new_kc;
}
void Virus::set_VtoD(double new_VtoD){
  _V_to_d = new_VtoD;
}





Virus::Virus(Virus* _parent, Host* _host, int _t, double _immK, double _tmpK){
  id = v_IDgenerator++;
  birth = _t;
  infectionK = _host->getInfectionHistory().size();
  bindingavid = bindingavid_ini = _parent->getBindingAvid();
  parent = _parent;
  level = parent->level + 1;
  distanceToParent = 0;
  if(_parent != NULL) distRoot = _parent->getDistRoot();
  else distRoot = 0;
  host = _host;
  immK = _immK;
  tmpK = _tmpK;
}

Virus::Virus(int _level, Virus* _parent, double _bindingavid, double _distance, Host* _host, int _t, double _immK, double _tmpK){
  id = v_IDgenerator++;
  birth = _t;
  infectionK = 0;
  parent = _parent;
  bindingavid = bindingavid_ini = _bindingavid;
  level = _level;
  distanceToParent = _distance;
  if(_parent != NULL) distRoot = _parent->getDistRoot() + _distance;
  else distRoot = 0;
  host = _host;
  immK = _immK;
  tmpK = _tmpK;
}


int Virus::getId(){
  return id;
}
double Virus::getImmK(){
  return immK;
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

int Virus::getHostK(){
  return(host->get_hostK());
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

Host* Virus::getHost(){
  return(host);
}

int Virus::getK(){
  return infectionK;
}

int Virus::getTmpImmK(){
  return tmpK;
}

void Virus::kill(int cur_t){
  death = cur_t;
}

double Virus::calculateRho(Host* _host){
  return(1 - pow((1/(_n*probSurvival(_host)*probReplication())),_v));
}

void Virus::mutate(){
  double bindingAvidChange;
  double change;
  double tmp = ((double) rand() / (RAND_MAX));
  std::exponential_distribution<double> dist(45);

  switch (_scenario){
  case 1: 
    if(tmp <= _prob_mut){
      change = dist(host->popn->generator);
      distanceToParent += change;
      distRoot += change;
    }
    break;
  case 2:
    bindingAvidChange = bindingavid_change(host);
    bindingavid += bindingAvidChange;
    break;
  case 3:
    if(tmp <= _prob_mut){
      change = dist(host->popn->generator);
      distanceToParent += change;
      distRoot += change;
    }
    bindingAvidChange = bindingavid_change(host);
    bindingavid += bindingAvidChange;
    distanceToParent += _V_to_d*fabs(bindingAvidChange);
    distRoot += _V_to_d*fabs(bindingAvidChange);
    break;
  case 4:
    if(tmp <= _prob_mut){
      change = dist(host->popn->generator);
      distanceToParent += change;
      distRoot += change;
    }
    bindingAvidChange = bindingavid_change(host);
    bindingavid += bindingAvidChange;
    break;
  default:
    break;
  }
}

double Virus::probSurvival(Host* _host){
  int size  = _host->getInfectionHistory().size();
  double tmp = 0;
  double tmp1 = 0;
  if(size > 0){
    tmp = getAntigenicDistance(this,_host->getInfectionHistory()[0]);
    for(int i = 1; i < size; ++i){
      tmp1 =getAntigenicDistance(this,_host->getInfectionHistory()[i]);
      if( tmp1 < tmp){
	tmp = tmp1;
      }
    }
  }
  tmpK = _r*_host->get_hostK() - 100*tmp;
  if(tmpK < 0){ tmpK = 0;}
  tmp = pow((1 - exp(-_p*(bindingavid + _q))), tmpK);
  return(tmp);
}
double Virus::getDistRoot(){
  return(distRoot);
}

double Virus::d_probSurvival(Host* _host){
  return(_p*immK*(pow((1-exp(-_p*(bindingavid+_q))),(immK-1)))*(exp(-_p*(bindingavid+_q))));
}

double Virus::probReplication(){
  return(exp(-_a*(pow(bindingavid,_b))));
}

double Virus::d_probReplication(){
  return((pow(-_a*_b*bindingavid,_b-1))*(exp(-_a*(pow(bindingavid,_b)))));
}	 

double Virus::bindingavid_change(Host* _host){
  double dV = probSurvival(_host)*d_probReplication() + probReplication()*d_probSurvival(_host);
  return(dV*_kc);
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
