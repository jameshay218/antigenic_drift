package ibid;

import java.util.ArrayList;
import java.util.List;

public class HostPopulation {
	private List<Host> susceptibles = new ArrayList<Host>();
	private List<Host> infecteds = new ArrayList<Host>();	
	private List<Host> recovereds = new ArrayList<Host>();		// this is the transcendental class, immune to all forms of virus  
	private int infectionsNow;
	private List<Host> newParents;
	private List<Host> newInfections;
	private List<Host> newRecoveries;
	private int day;
    public void stepForward(int d) {
    	day = d; 
    	//demographic
    	grow();
    	decline();
    	//--------------
    	infectionsNow = getI();
    	newParents = new ArrayList<Host>();    	// source hosts
    	newInfections = new ArrayList<Host>();	// challenging hosts
    	newRecoveries = new ArrayList<Host>();	// challenging hosts
    	contact();
    	recover();
    	wane();
    	updateStatus(); //update I and R status after contact and recovery is calculated 
    	mutateWithinHost();
    	//sample();
    	printStatus();
    }
	
    public HostPopulation() {
    	System.out.println("test");
    	for (int i = 0; i < Parameters.initialS; i++) {
			Host h = new Host();
			susceptibles.add(h);
		}
		
		// fill population with recovereds
		for (int i = 0; i < Parameters.initialR; i++) {
			Host h = new Host();			
			recovereds.add(h);
		}		
		
		for (int i = 0; i < Parameters.initialI; i++) {                    
			//Start with a homogenous viral population
			Virus v = new Virus();
            Host h = new Host(v);
            infecteds.add(h);	
		}
    }
 //
 	public int getN() {
 		return susceptibles.size() + infecteds.size() + recovereds.size();
 	}
 	public int getS() {
 		return susceptibles.size();
 	}
 	public int getI() {
 		return infecteds.size();
 	}
 	public int getR() {
 		return recovereds.size();
 	}	
 	public double getmeanTraitA() {
 		double meanTraitA = 0;
 		double tmptotalTraitA = 0;
 		for (int i=1;i<getI();i++) {
 			tmptotalTraitA = tmptotalTraitA + infecteds.get(i).getVirus().phenotype.getTraitA();
        }
 		meanTraitA = tmptotalTraitA/getI();
 		return meanTraitA;
 	}
 
 	public void grow() {
 		int newBirth = Random.nextPoisson(Parameters.mu*getN());
 		for (int i=1;i<newBirth;i++) { 
 			Host h = new Host();
 			susceptibles.add(h);
 		}
 	}
 	
 	public void decline() {
 		int newDeathS = Random.nextPoisson(Parameters.mu*getS());
 		for (int i=1;i<newDeathS;i++) {
 			susceptibles.remove(Random.nextInt(0, getS()-1));
 		}
 		int newDeathI = Random.nextPoisson(Parameters.mu*getI());
 		for (int i=1;i<newDeathI;i++) {
 			infecteds.remove(Random.nextInt(0, getI()-1));
 		}
 		int newDeathR = Random.nextPoisson(Parameters.mu*getR());
 		for (int i=1;i<newDeathR;i++) {
 			recovereds.remove(Random.nextInt(0, getR()-1));
 		}
 	}
 	
 	
 	public void contact() {


 	  int totalContacts = Random.nextPoisson(Parameters.contact*getI()*getS()/getN()); 
 	  //infective source
 	  for (int i=0;i<totalContacts;i++) {
 		  int sidx = Random.nextInt(0, getI()-1);
 		  newParents.add(infecteds.get(sidx)); //sampling with replacement
 	  }
 	  //System.out.println("pick infecteds source");
 	  
 	  for (int i=0;i<newParents.size();i++) {
 		 if(Random.nextBoolean(Parameters.rho)) {
 			 Virus v_src = newParents.get(i).getVirus();
 			 int hidx = Random.nextInt(0, getS()-1);
 			 Host h = susceptibles.get(hidx);  //sampling with replacement; java pass by address
 			 Virus v_tar = new Virus(v_src);
 			 h.addVirus(v_tar); //the virus has same trait as parent virus
 			 newInfections.add(h);
 			 //remove from susceptible
 			 susceptibles.remove(hidx);
 		 }
 	  }
 	  //infecteds.addAll(newInfections); //wait until
 	  //System.out.println("contact:" + totalContacts);
 	  //System.out.println("new I:" + newInfections.size());
 	}
 	
 	public void updateStatus() {
 		infecteds.addAll(newInfections);
 		recovereds.addAll(newRecoveries);
 	}
 	
 	public void mutateWithinHost() {
 	  for (int i=0;i<newInfections.size();i++) {
 		  //newInfections.get(i).getVirus().phenotype.getTraitA();
 		 newInfections.get(i).getVirus().phenotype.mutate();
 	  }
 	}
 	
 	public void recover() {
 		int totalRecovereds = Random.nextPoisson(Parameters.gamma*getI());
 		for (int i=0;i<totalRecovereds;i++) {
 		  int eidx = Random.nextInt(0, getI()-1);
 		  //recovereds.add(infecteds.get(eidx));
 		  newRecoveries.add(infecteds.get(eidx)); 
 		  infecteds.remove(eidx);
 		}
 		//System.out.println("recover:" + totalRecovereds);
 	}
 	
 	public void wane() {
 		int totalWanning = Random.nextPoisson(Parameters.wan*getR());
 		for (int i=0;i<totalWanning;i++) { //need t o write immune wanning to make SIR stay equilibrium 
 		  int widx = Random.nextInt(0, getR()-1);
 		  susceptibles.add(recovereds.get(widx));
 		  recovereds.remove(widx);
 	 	}
 	}
 	
 	public void printStatus() {
 	  System.out.println("Day\tS\tI\tR\tTraitA");
 	  System.out.println(day+"\t"+getS()+"\t"+getI()+"\t"+getR()+"\t"+getmeanTraitA());
 	}
}
