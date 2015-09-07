#include <iostream>
#include <vector>
#include "virus.hpp"
#include <time.h>

using namespace std;

int main(int argc, char *argv[]){
  vector<Virus> allViruses;
  Virus A(1, 0, NULL, 0.5, 0.5);
  allViruses.push_back(A);
  Virus B(2, &A, 0.5, 0.5);
  allViruses.push_back(B);
  Virus C(3, &A, 0.5, 0.3);
  allViruses.push_back(C);
  Virus D(4, &B, 0.5, 0.25);
  allViruses.push_back(D);
  Virus E(5, &B, 0.5, 0.4);
  allViruses.push_back(E);
  Virus F(6, &D, 0.5, 0.01);
  allViruses.push_back(F);

  Virus G(7, &C, 0.5, 2);
  allViruses.push_back(G);
  Virus H(8, &G, 0.5, 0.001);
  allViruses.push_back(H);
  Virus I(9, &G, 0.5, 0.002);
  allViruses.push_back(I);
  Virus J(10, &G, 0.5, 0.003);
  allViruses.push_back(J);

  cout << "------------- Antigenic Distance Test ----------------" << endl;
  cout << "Distance G to E: " << endl;
  cout << Virus::getAntigenicDistance(&G, &E) << endl;

  cout << "Distance F to A: " << endl;
  cout << Virus::getAntigenicDistance(&F, &A) << endl;

  cout << "Distance C to F: " << endl;
  cout << Virus::getAntigenicDistance(&C, &F) << endl;

  cout << "Distance C to E: " << endl;
  cout << Virus::getAntigenicDistance(&C, &E) << endl;
  
  cout << "Distance D to F: " << endl;
  clock_t t;
  t = clock();
  cout << "Time now: " << t << endl;
  for(int i = 0; i < 10000000; i++){
    Virus::getAntigenicDistance(&D, &F);
  }
  cout << "Time: " << clock() << endl;
  cout << "Time: " << (clock() - t)/CLOCKS_PER_SEC << endl;


  cout << "Distance H to J: " << endl;
  cout << Virus::getAntigenicDistance(&H, &J) << endl;

  cout << "Distance I to F: " << endl;
  cout << Virus::getAntigenicDistance(&I, &F) << endl;

  return 0;
}
