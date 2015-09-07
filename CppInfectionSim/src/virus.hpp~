#ifndef VIRUS_HPP
#define VIRUS_HPP

class Virus {
protected:
  int id;
  Virus* parent;
  double bindingavid;
  double distanceToParent;
  int level;
public:
  Virus();
  Virus(int _id, Virus* _parent, double _bindavid, double _distance);
  Virus(int _id, int _level, Virus* _parent, double _bindavid, double _distance);
  static double getAntigenicDistance(Virus* A, Virus* B);
  int getId();
  ~Virus(){};
  Virus* getParent(); 
  double getBindingAvid();
  double getDistance();
  int getLevel();
};

#endif
