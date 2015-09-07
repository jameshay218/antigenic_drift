#ifndef HOSTPOPULATION_HPP
#define HOSTPOPULATION_HPP

#include <vector>

class Host;
class VirusPopulation;

class HostPopulation{
private:
  std::vector<Host*> susceptibles;
  std::vector<Host*> infecteds;
  std::vector<Host*> recovereds;
  int day;
  double contactRate;

public:
  VirusPopulation* virusPopn;
  
  // Constructors
  HostPopulation();
  HostPopulation(int initialS, int initialI, int initialR, int iniDay, double _contactRate);
  ~HostPopulation(){};

  // Manage population temporal dynamics
  void stepForward();
  void grow();
  void decline();
  void contact();
  void recover();
  void wane();
  void setVirusPopn(VirusPopulation* _virusPopn);

  // Get properties of HostPopulation
  double getContactRate();
  int getDay();
  int countSusceptibles();
  int countInfecteds();
  int countRecovereds();
  
  // Print out current population status
  void printStatus();
  
};

#endif
