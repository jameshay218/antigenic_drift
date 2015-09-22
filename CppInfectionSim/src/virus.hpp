#ifndef VIRUS_HPP
#define VIRUS_HPP

class Host;

class Virus {
private:
  static int v_IDgenerator;

  int id;
  Virus* parent;
  Host* host;
  double bindingavid;
  double distanceToParent;
  int level;

public:
  static double _p;
  static double _r;
  static double _q;
  static double _a;
  static double _b;
  static double _n;
  static double _v;
  static double _k;
  
  static double getAntigenicDistance(Virus* A, Virus* B);

  // Constructors
  Virus();
  Virus(Virus* _parent, double _bindavid, double _distance, Host* _host);
  Virus(int _level, Virus* _parent, double _bindavid, double _distance, Host* _host);
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

};

#endif
