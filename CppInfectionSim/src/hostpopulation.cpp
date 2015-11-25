#include "hostpopulation.hpp"
#include "host.hpp"
#include "virus.hpp"

using namespace std;

HostPopulation::HostPopulation(){
  day = 1;
  contactRate = 0.5;
  mu = 1.0/(70.0*365.0);
  for(int i = 0; i < 70000; ++i){
    Host* firstH = new Host();
    susceptibles.push_back(firstH);
  }
    
}

HostPopulation::HostPopulation(int initialS, int initialI, int initialR, int iniDay, double _contactRate, double _mu, double _wane, double _gamma){
  Host* H;
  Virus* V;
  double iniBindingAvid = 0.5;

  day = iniDay;
  contactRate = _contactRate;
  mu = _mu;
  wane = _wane;
  gamma = _gamma;
  
  for(int i = 0; i < initialS;++i){
    H = new Host(Susceptible, this);
    susceptibles.push_back(H);
  }
  for(int i =0; i < initialI; ++i){
    H = new Host(Infected, this);
    V = new Virus(0, NULL, iniBindingAvid, 0, H,day);
    H->infect(V,day);
    infecteds.push_back(H);
  }
  for(int i = 0; i < initialR; ++i){
    H = new Host(Recovered, this);
    recovereds.push_back(H);
  }
}

HostPopulation::~HostPopulation(){
  cout << "Host Population Destructor" << endl;
  int j = susceptibles.size();
  for(int i = 0; i < j; ++i){
    delete susceptibles[i];
  }

  j = infecteds.size();
  for(int i = 0; i < j; ++i){
    delete infecteds[i];
  }

  j = recovereds.size();
  for(int i = 0; i < j; ++i){
    delete recovereds[i];
  }

  j = dead.size();
  for(int i = 0; i < j; ++i){
    delete dead[i];
  }
}


void HostPopulation::stepForward(int new_day){
  // Current day is changed
  day = new_day;

  //New births of susceptibles
  grow();

  // Deaths from all compartments
   decline();

  // Infected population grows from transmission
  contact();

  // Infected population declines from recovery
  recoveries();

  // Some immunity wanes
   waning();
  
   // Mutations?
   mutations();
   // Go through each host in infected vector and mutate its virus

}


void HostPopulation::grow(){
  poisson_distribution<int> poisson(mu*countN());
  int newBirths = poisson(generator);
  for(int i = 0; i < newBirths; ++i){
    Host* h = new Host(Susceptible, this);
    susceptibles.push_back(h);
  }
}

void HostPopulation::decline(){
  poisson_distribution<int> poissonS(mu*countSusceptibles());
  poisson_distribution<int> poissonI(mu*countInfecteds());
  poisson_distribution<int> poissonR(mu*countRecovereds());
  int newDeaths = poissonS(generator);
  int index;
  for(int i = 0; i < newDeaths; ++i){
    index = rand() % countSusceptibles();
    dead.push_back(susceptibles[index]);
    susceptibles[index]->die(day);
    susceptibles.erase(susceptibles.begin()+index);
  }

  newDeaths = poissonI(generator);
  for(int i = 0; i < newDeaths; ++i){
    index = rand() % countInfecteds();
    dead.push_back(infecteds[index]);
    infecteds[index]->die(day);
    infecteds.erase(infecteds.begin()+index);
  }

  newDeaths = poissonR(generator);
  for(int i = 0; i < newDeaths; ++i){
    index = rand() % countRecovereds();
    dead.push_back(recovereds[index]);
    recovereds[index]->die(day);
    recovereds.erase(recovereds.begin()+index);
  }
}

void HostPopulation::contact(){
  poisson_distribution<int> poisson(contactRate*countInfecteds()*countSusceptibles()/countN());
  int totalContacts = poisson(generator);
  int index1 = 0;
  int index2 = 0;
  double tmp = 0;
  for(int i = 0; i < totalContacts; ++i){
    index1 = rand() % countInfecteds();
    tmp = ((double) rand() / (RAND_MAX));
    index2 = rand() % countSusceptibles();
    //cout << "Prob survival: " << infecteds[index1]->getCurrentVirus()->calculateRho(susceptibles[i]) << endl;
    if(tmp <= infecteds[index1]->getCurrentVirus()->calculateRho(susceptibles[index2])){
      Virus* newV = new Virus(infecteds[index1]->getCurrentVirus(), susceptibles[index2], day);
      susceptibles[index2]->infect(newV,day);
      infecteds.push_back(susceptibles[index2]);
      susceptibles.erase(susceptibles.begin() + index2);
    }
  }
}

void HostPopulation::recoveries(){
  poisson_distribution<int> poisson(gamma*countInfecteds());
  int noRecovered = poisson(generator);
  int index = 0;
  for(int i = 0; i < noRecovered; ++i){
    index = rand() % countInfecteds();
    recovereds.push_back(infecteds[index]);
    infecteds[index]->recover(day);
    infecteds.erase(infecteds.begin()+index);    
  }
}

void HostPopulation::waning(){
  poisson_distribution<int> poisson(wane*countRecovereds());
  int noWane = poisson(generator);
  int index = 0;
    for(int i = 0; i < noWane; ++i){
    index = rand() % countRecovereds();
    susceptibles.push_back(recovereds[index]);
    recovereds[index]->wane();
    recovereds.erase(recovereds.begin()+index);    
    }
}

void HostPopulation::mutations(){
  int j = infecteds.size();
  for(int i = 0; i < j; ++i){
    infecteds[i]->getCurrentVirus()->mutate();
  }
}

int HostPopulation::countSusceptibles(){
  return(susceptibles.size());
}

int HostPopulation::countInfecteds(){
  return(infecteds.size());
}

int HostPopulation::countRecovereds(){
  return(recovereds.size());
}

int HostPopulation::getDay(){
  return(day);
}

int HostPopulation::countN(){
  return(infecteds.size() + recovereds.size() + susceptibles.size());
}

double HostPopulation::getContactRate(){
  return(contactRate);
}


void HostPopulation::printStatus(){
  cout << "Current day: " << day << endl;
  cout << "Susceptible: " << susceptibles.size() << endl;
  cout << "Infecteds: " << infecteds.size() << endl;
  cout << "Recovereds: " << recovereds.size() << endl;
  cout << "Died: " << dead.size() << endl;
  cout << "Total: " << countN() << endl;
}


void HostPopulation::writeViruses(std::ofstream& output, std::string filename){
  output.open(filename);
  int x = 0;
  vector<vector<Host*>> tmp = {susceptibles, infecteds, recovereds, dead};
  int q = tmp.size();  
  int j = 0;
  output << "vid, birth, death, parent, infectionK, bindingAvid, distance" << endl;
  for(int y = 0; y < q; ++y){
    j = tmp[y].size();
    for(int i = 0; i < j; ++i){
      x = tmp[y][i]->getInfectionHistory().size();
      if(x > 0){
	for(int ii = 0; ii < x; ++ii){
	  output << tmp[y][i]->getInfectionHistory()[ii]->getId() << ",";
	  output << tmp[y][i]->getInfectionHistory()[ii]->getBirth() << ",";
	  output << tmp[y][i]->getInfectionHistory()[ii]->getDeath() << ",";
	  if(tmp[y][i]->getInfectionHistory()[ii]->getParent()){
	    output << tmp[y][i]->getInfectionHistory()[ii]->getParent()->getId() << ",";
	  } else {
	    output << 0 << ",";
	  }
	  output << tmp[y][i]->getInfectionHistory()[ii]->getK() << ",";
	  output << tmp[y][i]->getInfectionHistory()[ii]->getBindingAvid() << ",";
	  output << tmp[y][i]->getInfectionHistory()[ii]->getDistance() << endl;
	}
      }
    }
  }
  output.close();
}
