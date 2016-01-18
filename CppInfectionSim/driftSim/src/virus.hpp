#ifndef VIRUS_HPP
#define VIRUS_HPP

class Host;

class Virus {
private:
  static int v_IDgenerator;

  int id;
  int birth;
  int death;
  int infectionK;
  Virus* parent;
  Host* host;
  double bindingavid_ini;
  double bindingavid;
  double distanceToParent;
  double immK;
  double tmpK;
  double distRoot;
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
  static int _scenario;

  static double getAntigenicDistance(Virus* A, Virus* B);

  // Constructors
  Virus();
  Virus(Virus* _parent, Host* _host, int _t, double _immK, double _tmpK);
  Virus(int _level, Virus* _parent, double _bindingavid, double _distance, Host* _host, int _t, double _immK, double _tmpK);
  ~Virus(){};  // Don't really need to worry about pointers. VirusPopulation should take care of memory.

  // Accessing attributes
  int getId();
  int getBirth();
  int getDeath();
  int getK();
  int getTmpImmK();
  int getHostK();
  Host* getHost();

  double getImmK();
  Virus* getParent(); 
  double getBindingAvid();
  double getIniBindingAvid();
  double getDistance();
  double getDistRoot();
  int getLevel();

  // Calculations/events
  double calculateRho(Host* _host);
  void mutate();
  double probSurvival(Host* _host);
  double probReplication();
  double d_probSurvival(Host* _host);
  double d_probReplication();
  double bindingavid_change(Host* _host);
  void kill(int cur_t);

  // Change static member variables
  static void set_default();
  static void set_p(double _new_p);
  static void set_r(double new_r);
  static void set_q(double new_q);
  static void set_a(double new_a);
  static void set_b(double new_b);
  static void set_n(double new_n);
  static void set_v(double new_v);
  static void set_prob_mut(double new_probMut);
  static void set_exp_dist(double new_exp);
  static void set_kc(double new_kc);
  static void set_VtoD(double new_VtoD);
  static void set_scenario(int _scen);
};

#endif
