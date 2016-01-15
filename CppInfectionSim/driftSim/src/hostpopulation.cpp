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

HostPopulation::HostPopulation(int initialS, int initialI, int initialR, int iniDay, double _contactRate, double _mu, double _wane, double _gamma, double _iniBindingAvid){
  Host* H;
  Virus* V;
  double iniBindingAvid = _iniBindingAvid;
  int _k = 0;
  generator.seed(time(NULL));
  day = iniDay;
  contactRate = _contactRate;
  mu = _mu;
  wane = _wane;
  gamma = _gamma;
  
  for(int i = 0; i < initialS;++i){
    //  _k = rand() % 20;
    H = new Host(Susceptible, this, _k);
    susceptibles.push_back(H);
  }
  for(int i =0; i < initialI; ++i){
    H = new Host(Infected, this);
    V = new Virus(0, NULL, iniBindingAvid, 0, H,day,0,0);
    H->infect(V,day);
    infecteds.push_back(H);
  }
  for(int i = 0; i < initialR; ++i){
    H = new Host(Recovered, this);
    recovereds.push_back(H);
  }
}

HostPopulation::~HostPopulation(){
  // Rcpp::Rcout << "Host Population Destructor" << endl;
  int j = susceptibles.size();
  for(int i = 0; i < j; ++i){
    //Rcpp::Rcout << i << endl;
    //    Rcpp::Rcout << susceptibles[i]->getState() << endl;
    delete susceptibles[i];
  }
  //Rcpp::Rcout << "S deleted" << endl;
  j = infecteds.size();
  for(int i = 0; i < j; ++i){
    delete infecteds[i];
  }
  // Rcpp::Rcout << "I deleted" << endl;
  j = recovereds.size();
  for(int i = 0; i < j; ++i){
    delete recovereds[i];
  }
  //  Rcpp::Rcout << "R deleted" << endl;
  j = dead.size();
  for(int i = 0; i < j; ++i){
    // if(!dead[i]->isDead()){
      // Rcpp::Rcout << dead[i]->getState() << endl;
    //  }
    delete dead[i];
  }
  // Rcpp::Rcout << "D deleted" << endl;
}

void HostPopulation::stepForward(int new_day){
  // Current day is changed
  day = new_day;

  //New births of susceptibles
  // Rcpp::Rcout << "Grow" << endl;
  grow();

  // Infected population grows from transmission
  //Rcpp::Rcout << "Contact" << endl;
  contact();

  // Infected population declines from recovery
  // Rcpp::Rcout << "Recovered" << endl;
  recoveries();

  // Some immunity wanes
  // Rcpp::Rcout << "Waning" << endl;
  waning();

  // Deaths from all compartments - MUST BE LAST EVENT
  //Rcpp::Rcout << "Decline" << endl;
  decline();
  
  //Rcpp::Rcout << "Update" << endl;
  updateCompartments();

  // Go through each host in infected vector and mutate its virus
  // Mutations?
  mutations();
  

}


void HostPopulation::grow(){
  /*  poisson_distribution<int> poisson(mu*countN());
      int newBirths = poisson(generator);*/
  int newBirths = R::rpois(mu*countN());
  for(int i = 0; i < newBirths; ++i){
    Host* h = new Host(Susceptible, this);
    new_births.push_back(h);
  }
}

void HostPopulation::decline(){
  /*poisson_distribution<int> poissonS(mu*countSusceptibles());
  poisson_distribution<int> poissonI(mu*countInfecteds());
  poisson_distribution<int> poissonR(mu*countRecovereds());
  int newDeaths = poissonS(generator);*/
  int newDeaths = R::rpois(mu*countSusceptibles());
  int totDeaths = 0;
  int index;
  for(int i = 0; i < newDeaths; ++i){
    if(countSusceptibles() > 0){
      index = floor(R::unif_rand()*(countSusceptibles()));
      dead.push_back(susceptibles[index]);
      susceptibles[index]->die(day);
      susceptibles[index] = susceptibles.back();
      susceptibles.pop_back();
      totDeaths++;
    }
  }
  newDeaths = R::rpois(mu*countInfecteds());
  //  newDeaths = poissonI(generator);
  for(int i = 0; i < newDeaths; ++i){
    if(countInfecteds() > 0){
      index = floor(R::unif_rand()*(countInfecteds()));
      dead.push_back(infecteds[index]);
      infecteds[index]->die(day);
      infecteds[index] = infecteds.back();
      infecteds.pop_back();
      totDeaths++;
    }
  }
  newDeaths = R::rpois(mu*countRecovereds());
  //  newDeaths = poissonR(generator);
  for(int i = 0; i < newDeaths; ++i){
    if(countRecovereds() > 0){
      index = floor(R::unif_rand()*(countRecovereds()));
      dead.push_back(recovereds[index]);
      recovereds[index]->die(day);
      recovereds[index] = recovereds.back();
      recovereds.pop_back();
      totDeaths++;
    }
  }
}

void HostPopulation::contact(){
  // Generate number of contacts between infecteds and susceptibles
  //poisson_distribution<int> poisson(contactRate*countInfecteds()*countSusceptibles()/countN());
  //  int totalContacts = poisson(generator);
  int totalContacts = R::rpois(contactRate*countInfecteds()*countSusceptibles()/countN());
  int index1 = 0;
  int index2 = 0;
  double tmp = 0;
  double number_success = 0;
  
  /* For each contact, get a random I and random S. With a probability proportional to virus survival, infect the susceptible. 
     Note that the susceptible may contact multiple infecteds */
  for(int i = 0; i < totalContacts; ++i){
    if(countInfecteds() > 0){
      index1 = floor(R::unif_rand()*(countInfecteds()));
      index2 = floor(R::unif_rand()*(countSusceptibles()));
      //      Rcpp::Rcout << index1 << " and " << index2 << endl;
      //      tmp = ((double) rand() / (RAND_MAX));
      tmp = R::unif_rand();
      // Rcpp::Rcout << "Prob survival: " << infecteds[index1]->getCurrentVirus()->calculateRho(susceptibles[i]) << endl;
      if(tmp <= infecteds[index1]->getCurrentVirus()->calculateRho(susceptibles[index2])){
	if(susceptibles[index2]->isSusceptible()){
	  Virus* newV = new Virus(infecteds[index1]->getCurrentVirus(), susceptibles[index2], day, infecteds[index1]->getCurrentVirus()->getTmpImmK(),susceptibles[index2]->getInfectionHistory().size()+1);
	  susceptibles[index2]->infect(newV,day);
	  new_infecteds.push_back(susceptibles[index2]); // Record new infected to add later
	  number_success++;
	}
      }
    }
  }
  //  Rcpp::Rcout << "Beta: " << number_success/totalContacts << endl;
}

void HostPopulation::recoveries(){
  // Random number of recoveries of current infecteds
  /*  poisson_distribution<int> poisson(gamma*countInfecteds());
      int noRecovered = poisson(generator);*/
  int noRecovered = R::rpois(gamma*countInfecteds());
  int index = 0;
  // For each recovery, choose a random infected individual and add them to the list of newly recovered individuals
  for(int i = 0; i < noRecovered; ++i){
    if(countInfecteds() > 0){
      index = floor(R::unif_rand()*(countInfecteds()));
      // Record the new recovered individual and temporarily delete from the list of infecteds
      new_recovereds.push_back(infecteds[index]);
      infecteds[index]->recover(day);
      infecteds[index] = infecteds.back();
      infecteds.pop_back();
    }
  }
  // Add the newly recovered individuals back to list of infecteds so that other events can be calculated correctly.
  infecteds.insert(infecteds.end(),new_recovereds.begin(),new_recovereds.end());
}

void HostPopulation::waning(){
  /*  poisson_distribution<int> poisson(wane*countRecovereds());
      int noWane = poisson(generator);*/
  int noWane = R::rpois(wane*countRecovereds());
  int index = 0;
  for(int i = 0; i < noWane; ++i){
    if(countRecovereds() > 0){
      index = floor(R::unif_rand()*(countRecovereds()));
      new_susceptibles.push_back(recovereds[index]);
      recovereds[index]->wane();
      recovereds[index] = recovereds.back();
      recovereds.pop_back();
    }
  }
  // Add the newly susceptible individuals back to list of recovereds so that other events can be calculated correctly.
  recovereds.insert(recovereds.end(),new_susceptibles.begin(),new_susceptibles.end());
}

void HostPopulation::updateCompartments(){
  infecteds.insert(infecteds.end(),new_infecteds.begin(),new_infecteds.end());
  recovereds.insert(recovereds.end(),new_recovereds.begin(),new_recovereds.end());
  susceptibles.insert(susceptibles.end(),new_susceptibles.begin(),new_susceptibles.end());
  susceptibles.insert(susceptibles.end(),new_births.begin(),new_births.end());

  susceptibles.erase(remove_if(susceptibles.begin(),susceptibles.end(), [](Host *h){return(!h->isSusceptible());}),susceptibles.end());
  infecteds.erase(remove_if(infecteds.begin(),infecteds.end(), [](Host* o){return(!o->isInfected());}),infecteds.end());
  recovereds.erase(remove_if(recovereds.begin(),recovereds.end(), [](Host* o){return(!o->isRecovered());}),recovereds.end());

  new_infecteds.clear();
  new_recovereds.clear();
  new_susceptibles.clear();
  new_births.clear();
}

void HostPopulation::mutations(){
  int j = infecteds.size();
  if(j > 0){
    for(int i = 0; i < j; ++i){
      infecteds[i]->getCurrentVirus()->mutate();
    }
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
  Rcpp::Rcout << "Current day: " << day << endl;
  Rcpp::Rcout << "Susceptible: " << susceptibles.size() << endl;
  Rcpp::Rcout << "Infecteds: " << infecteds.size() << endl;
  Rcpp::Rcout << "Recovereds: " << recovereds.size() << endl;
  Rcpp::Rcout << "Died: " << dead.size() << endl;
  Rcpp::Rcout << "Total: " << countN() << endl;
}


void HostPopulation::writeViruses(std::ofstream& output, std::string filename){
  Rcpp::Rcout << "#########################" << endl;
  Rcpp::Rcout << "Writing viruses to csv..." << endl;
  output.open(filename);
  int x = 0;
  vector<vector<Host*>> tmp = {susceptibles, infecteds, recovereds, dead};
  vector<Virus*> viruses;
  int q = tmp.size();  
  int j = 0;
  for(int y = 0; y < q; ++y){
    j = tmp[y].size();
    for(int i = 0; i < j; ++i){
      x = tmp[y][i]->getInfectionHistory().size();
      if(x > 0){
	for(int ii = 0; ii < x; ++ii){
	  viruses.push_back(tmp[y][i]->getInfectionHistory()[ii]);
	}
      }
    }
  }
  Rcpp::Rcout << "Sorting viruses..." << endl;
  sort(viruses.begin(),viruses.end(),[](Virus* lhs, Virus* rhs){
      return(lhs->getId() < rhs->getId());
    });
  Rcpp::Rcout << "Writing output..." << endl;
  output << "vid, birth, death, parentid, bindingAvidityIni, bindingAvidityFinal, infection_no, distance_to_parent, hostK, host_infections, immK, distRoot" << endl;
  //  output << "hostK, host_infections, vid, birth, death, parentid, infection_no, bindingAvidIni, bindingAvidFinal, distance_to_parent, immK, distRoot" << endl;
  j = viruses.size();
  for(int i = 0; i < j;++i){
    output << viruses[i]->getId() << ",";
    output << viruses[i]->getBirth() << ",";
    output << viruses[i]->getDeath() << ",";
    if(viruses[i]->getParent()){
      output << viruses[i]->getParent()->getId() << ",";
    } else {
      output << 0 << ",";
    }
    output << viruses[i]->getIniBindingAvid() << ",";
    output << viruses[i]->getBindingAvid() << ",";
    output << viruses[i]->getK() << ",";
    output << viruses[i]->getDistance() << ",";
    output << viruses[i]->getHostK() << ",";
    output << viruses[i]->getHost()->getInfectionHistory().size() << ",";
    output << viruses[i]->getImmK() << ",";
    output << viruses[i]->getDistRoot() << endl;
  }

  output.close();
  Rcpp::Rcout << "Writing viruses complete" << endl;
  Rcpp::Rcout << "#########################" << endl;
  Rcpp::Rcout << endl;
}

void HostPopulation::virusPairwiseMatrix(std::ofstream& output, std::string filename, int sampSize){
  Rcpp::Rcout << "#########################" << endl;
  Rcpp::Rcout << "Writing pairwise matrix" << endl;
  vector<Host*> hosts;
  vector<Virus*> viruses;
  vector<int> sampIndices;
  output.open(filename);
  hosts.insert(hosts.end(),susceptibles.begin(),susceptibles.end());
  hosts.insert(hosts.end(),infecteds.begin(),infecteds.end());
  hosts.insert(hosts.end(),recovereds.begin(),recovereds.end());
  hosts.insert(hosts.end(),dead.begin(),dead.end());
  int j = hosts.size();
  int x = 0;
  Rcpp::Rcout << "Accumulating all viruses..." << endl;
  for(int i = 0; i < j; ++i){
    x = hosts[i]->getInfectionHistory().size();
    if(x > 0){
      for(int ii = 0; ii < x; ++ii) {
	viruses.push_back(hosts[i]->getInfectionHistory()[ii]);
      }
    }
  }
  Rcpp::Rcout << "Total number of viruses: " << viruses.size() << endl;
  j = viruses.size();
  for(int i = 0; i < j;++i){
    sampIndices.push_back(i);
  }
  random_shuffle(sampIndices.begin(),sampIndices.end());
  sampIndices.resize(sampSize);
  sort(sampIndices.begin(),sampIndices.end());  

  Rcpp::Rcout << "Sorting virus vector..." << endl;

  sort(viruses.begin(),viruses.end(),[](Virus* lhs, Virus* rhs){
      return(lhs->getId() < rhs->getId());
    });
  
  Rcpp::Rcout << "Finding pairwise matrix..." << endl;
 
  j = sampIndices.size();
  int index = 0;
  int index2 = 0;
  for(int row = 0; row < j; ++row){
    index = sampIndices[row];
    for(int col = 0; col < row-1; ++col){
      output << ",";
    }
    for(int col = row; col < j; ++col){
      output << ",";
      index2 = sampIndices[col];
      output << Virus::getAntigenicDistance(viruses[index],viruses[index2]);
    }
    output << endl;
  }
  Rcpp::Rcout << "Pairwise output complete" << endl;
  Rcpp::Rcout << "#########################" << endl;
}  
  
