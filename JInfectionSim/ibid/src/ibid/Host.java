package ibid;
/* A human individual that harbors viruses and immunity */

import java.util.*;
import java.io.*;
import java.util.regex.*;

public class Host {

	// fields
	//private Virus expose;    // save current exposed virus  
	private Virus infection; // save current infected virus												
	private Phenotype[] immuneHistory = new Phenotype[0];
	//private int isExposed;
	private int isInfected;
	
	
	// naive host
	public Host() {
		//initializeHistory();		
	}
	
	// initial infected host
	public Host(Virus v) {
		//expose = v;
		infection = v;
		isInfected = 1;
		//isExposed = 1;
		//initializeHistory();
	}
	
	
	public void onSet() {
		//should check risk first
		//infection = expose;
		//isExposed = 0;
		//isInfected = 1;
	}
	
	public Virus getVirus() {
		return infection;
	}
	
	public void addVirus(Virus v) {
		infection = v;
		isInfected = 1;
	}
	
	// cross immunity between a virus phenotype and a host's immune history
	// here encoded more directly as risk of infection, which ranges from 0 to 1
	public double getRiskOfInfection( Phenotype[] history) {
	
		// find closest phenotype in history
		double closestDistance = 100.0;
		if (history.length > 0) {
			for (int i = 0; i < history.length; i++) {
				double thisDistance = infection.phenotype.getTraitDistance(history[i]);
				if (thisDistance < closestDistance) {
					closestDistance = thisDistance;
				}
				if (thisDistance < 0.01) {
					break;
				}
			}
		} 
		
		double risk = closestDistance;
		return risk;
		
	}
}
