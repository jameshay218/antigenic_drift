#ifndef HOST_HPP
#define HOST_HPP

#include <vector>

enum State {Susceptible, Infected, Recovered};

class Virus;

class Host{
private:
  State state;
  std::vector<Virus*> infectionHistory;
  Virus* currentInfection;

public:
  // Constructors
  Host();
  Host(State _state);
  ~Host(){};

  // Calculations/Events
  double calculateBeta();
  void infect(Virus* newInfection);
  void recover();

  // Attribute Access
  Virus* getCurrentVirus();
  std::vector<Virus*> getInfectionHistory();
  State getState();

};

#endif
