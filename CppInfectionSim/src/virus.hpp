#ifndef VIRUS_HPP
#define VIRUS_HPP

class Host;

class Virus {
private:
  static int v_IDgenerator;

  int id;
  int birth;
  int death = -1;
  int infectionK;
  Virus* parent;
  Host* host;
  double bindingavid_ini;
  double bindingavid;
  double distanceToParent;
  double immK;
  double tmpK;
  int level;



public:
  static double _p;
  static double _r;
  static double _q;
  static double _a;
  static double _b;
  static double _n;
  static double _v;
  static double _prob_mut;
  static double _exp_dist;
  static double _kc;
  static double _V_to_d;

  static double getAntigenicDistance(Virus* A, Virus* B);

  // Constructors
  Virus();
  Virus(Virus* _parent, Host* _host, int _t);
  Virus(int _level, Virus* _parent, double _bindingavid, double _distance, Host* _host, int _t);
  ~Virus(){};  // Don't really need to worry about pointers. VirusPopulation should take care of memory.

  // Accessing attributes
  int getId();
  int getBirth();
  int getDeath();
  int getK();
  double getImmK();
  Virus* getParent(); 
  double getBindingAvid();
  double getIniBindingAvid();
  double getDistance();
  int getLevel();

  // Calculations/events
  void updateK(int _k);
  double calculateRho(Host* _host);
  void mutate();
  double probSurvival(Host* _host);
  double probReplication();
  double d_probSurvival(Host* _host);
  double d_probReplication();
  double bindingavid_change(Host* _host);
  void kill(int cur_t);

};

#endif
