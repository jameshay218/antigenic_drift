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

  double hostK;

public:
  HostPopulation* popn;
  // Constructors
  Host();
  Host(State _state, HostPopulation* _popn);
  Host(State _state, HostPopulation* _popn, int _k);
  ~Host();

  // Calculations/Events
  double calculateBeta();
  void infect(Virus* newInfection, int cur_t);
  void recover(int cur_t);
  void wane();
  void die(int cur_t);
  double get_hostK();

  bool isInfected();
  bool isDead();
  bool isSusceptible();
  bool isRecovered();

  // Attribute Access
  Virus* getCurrentVirus();
  std::vector<Virus*> getInfectionHistory();
  State getState();
  std::default_random_engine get_generator();

};

#endif
