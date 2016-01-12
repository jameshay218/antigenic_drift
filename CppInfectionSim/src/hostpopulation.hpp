#ifndef HOSTPOPULATION_HPP
#define HOSTPOPULATION_HPP

#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

class Host;

class HostPopulation{
private:
  std::vector<Host*> susceptibles;
  std::vector<Host*> infecteds;
  std::vector<Host*> recovereds;
  std::vector<Host*> dead;

  std::vector<Host*> new_infecteds;
  std::vector<Host*> new_recovereds;
  std::vector<Host*> new_susceptibles;
  std::vector<Host*> new_births;

  int day;
  double contactRate;
  double mu;
  double wane;
  double gamma;


public:
  std::default_random_engine generator;

  // Constructors
  HostPopulation();
  HostPopulation(int initialS, int initialI, int initialR, int iniDay, double _contactRate, double _mu, double _wane, double _gamma, double _iniBindingAvid);
  ~HostPopulation();

  // Manage population temporal dynamics
  void stepForward(int new_day);
  void grow();
  void decline();
  void contact();
  void recoveries();
  void waning();
  void mutations();
  void updateCompartments();

  // Get properties of HostPopulation
  double getContactRate();
  int getDay();
  int countSusceptibles();
  int countInfecteds();
  int countRecovereds();
  int countN();
  
  // Print out current population status
  void printStatus();
  void writeViruses(std::ofstream& output, std::string filename);
  void virusPairwiseMatrix(std::ofstream& output, std::string filename, int sampSize);
};

#endif
