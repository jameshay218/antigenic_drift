#ifndef HOST_HPP
#define HOST_HPP

#include <vector>
#include<cstdlib>

#include "hostpopulation.hpp"
#include "virus.hpp"

enum State {Susceptible, Infected, Recovered, Dead};

class Host{
private:
  State state;
  std::vector<Virus*> infectionHistory;
  Virus* currentInfection;
  HostPopulation* popn;

public:
  // Constructors
  Host();
  Host(State _state, HostPopulation* _popn);
  ~Host();

  // Calculations/Events
  double calculateBeta();
  void infect(Virus* newInfection, int cur_t);
  void recover(int cur_t);
  void wane();
  void die(int cur_t);

  // Attribute Access
  Virus* getCurrentVirus();
  std::vector<Virus*> getInfectionHistory();
  State getState();

};

#endif
