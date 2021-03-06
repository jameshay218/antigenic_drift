package ibid;

import java.util.*;

public class Simulation {
    public List<Integer> list1;
    public List<Integer> list2; 
	
    public void test() {
    	//for (int i=0;i<200;i++) {
    	//System.out.print(""+Random.nextInt(0,1));
    	//}
    	
    	 double totalContacts = 0;
    	 for (int i=0;i<10000;i++) {
    	    totalContacts = totalContacts + Random.nextPoisson(500);
    	  }
    	  double meanContacts = totalContacts/10000;
    	  System.out.println(meanContacts);
    }
    public void run() {
    	//HostPopulation hpo = new HostPopulation();
    	int d = 1;
    	HostPopulation hpop = new HostPopulation();
    	hpop.printStatus();
    	while(d<Parameters.day) {
    		hpop.stepForward(d);
    	    d++;
    	}
    	  //int passing = 3;

    	  //Receiving (passing);
    	 
    	  //System.out.println("The value of passing is: " + passing);
	    System.out.println(list1); //prints 1, 2, 3
	    swapList(list1, list2);
	    System.out.println(list1); //prints 1, 2, 3
	    System.out.println(list2); 
    }
    
    Simulation() {
	    list1 = new ArrayList<Integer>(Arrays.asList(1, 2, 3));
	    list1.add(3);
	    list2 = new ArrayList<Integer>(Arrays.asList(4, 5, 6));
    }
    
    public static void main(String[] args) {
		
    	Simulation s = new Simulation();
        s.run();

		
		//Random Random()  
		//Creates a new random number generator.
		//cern.jet.random.AbstractDistribution.makeDefaultGenerator();
		// TODO Auto-generated method stub
		//HostPopulation hp = new HostPopulation();
		//for (int i=0;i<10;i++) {
		//	hp.stepForward();
    	//}
		
		//double mu = Parameters.getBirth();
		//double totalContactRate = 1;
        //int contacts = Random.nextPoisson(totalContactRate);
        //System.out.println(contacts);
	}
	public void swapList(List<Integer> l1, List<Integer> l2){
	    List<Integer> tmpList = l1;
	    tmpList.add(2);
	    tmpList.add(1);
	    this.list1 = this.list2;
	    this.list2 = tmpList;
	    //tmpList.clear();
	}
	
	public void Receiving (int var)
	{
	  var = var + 2;
	}
	
}
