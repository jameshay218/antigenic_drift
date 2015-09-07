#ifndef VIRUSPOPULATION_HPP
#define VIRUSPOPULATION_HPP

#include <vector>
#include <cstddef>
#include "virus.hpp"

class VirusPopulation{
private:
double p;
double r;
double q;
double a;
double b;
double n;
std::vector<Virus*> viruses;

public:
// Constructors
VirusPopulation();
VirusPopulation(double _p, double _r, double _q, double _a, double _b, double _n);
~VirusPopulation(){};

// Calculations/Events
void addVirus(Virus* newVirus);
double getAntigenicDistance(Virus* A, Virus* B);

// Attribute access
double getP(){return p;};
double getR(){return r;};
double getQ(){return q;};
double getA(){return b;};
double getB(){return b;};
double getN(){return n;};

};

#endif
