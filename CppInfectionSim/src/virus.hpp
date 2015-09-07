#ifndef VIRUS_HPP
#define VIRUS_HPP

class Host;
class VirusPopulation;

class Virus {
protected:
  int id;
  Virus* parent;
  Host* host;
  VirusPopulation* population;
  double bindingavid;
  double distanceToParent;
  int level;

public:
  // Constructors
  Virus();
  Virus(int _id, Virus* _parent, double _bindavid, double _distance, Host* _host, VirusPopulation* popn);
  Virus(int _id, int _level, Virus* _parent, double _bindavid, double _distance, Host* _host, VirusPopulation* popn);
  ~Virus(){};  // Don't really need to worry about pointers. VirusPopulation should take care of memory.

  // Accessing attributes
  int getId();
  Virus* getParent(); 
  double getBindingAvid();
  double getDistance();
  int getLevel();

  // Calculations/events
  double calculateRho();
  void mutate();
  double probSurvival();
  double probReplication();

  // Destructor


};

#endif
