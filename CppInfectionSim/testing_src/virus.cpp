#include "virus.hpp"
#include<iostream>
using namespace std;

Virus::Virus(int _id, Virus* _parent, double _bindavid, double _distance){
  id = _id;
  bindingavid = _bindavid;
  parent = _parent;
  level = parent->level + 1;
  distanceToParent = _distance;
}


Virus::Virus(int _id, int _level, Virus* _parent, double _bindavid, double _distance){
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
}
