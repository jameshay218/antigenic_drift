package ibid;
/* Antigenic phenotype present in individual Viruses and within Hosts as immune history */
/* Should be able to calculate distance and cross-immunity between two phenotypes */
/* Moving up to multiple dimensions is non-trivial and requires thought on the implementation */
/* Multiple Viruses can reference a single Phenotype object */

import static java.lang.Math.*;

import java.util.*;

public class Phenotype {

	// fields
	private double traitA;
	private String traitAName;
	private double traitB;
	private String traitBName;
	private double traitC;
	private String traitCName;
	
	// constructor
	public Phenotype() {
	  traitA = 0; //reset to 0 if there are no antigenic change
	  traitAName = "Antigenic Change"; 
	  traitB = 0;
	  traitBName = "undefined";
	  traitC = 0;
	  traitCName = "undefined";
	}
	// constructor
	public Phenotype(Phenotype another) {
	  //traitA = another.traitA;
	  traitA = 0;
	  traitAName = another.traitAName;
	  traitB = another.traitB;
	  traitBName = another.traitBName;
	  traitC = another.traitC;
	  traitCName = another.traitCName;
	}
	
	public Phenotype(double tA, double tB, double tC) {
		traitA = tA;
		traitB = tB;
		traitC = tC;
	}
		
	public double getTraitA() {
		return traitA;
	}
	public double getTraitB() {
		return traitB;
	}
	public double getTraitC() {
		return traitC;
	}	
	
	public void setTraitA(double tA) {
		traitA = tA;
	}
	public void setTraitB(double tB) {
		traitB = tB;
	}
	public void setTraitC(double tC) {
		traitC = tC;
	}	
		
	// raw antigenic distance between two phenotypes
	public double getTraitDistance(Phenotype p) {
		Phenotype p2d = p;
		double distA = (getTraitA() - p2d.getTraitA());
		double distB = (getTraitB() - p2d.getTraitB());	
		double dist = (distA * distA) + (distB * distB);
		dist = Math.sqrt(dist);
		System.out.println(dist);
		return dist;
	}


	
	
	// returns a mutated copy, original Phenotype is unharmed
	// mutate with gamma
	//public Phenotype mutate() {
	//		Phenotype p = new Phenotype();
	//		return p;
	//}
	
	public void mutate() {
		traitA = 1;
		System.out.println(traitA);
    }	

	
	
	public String toString() {
		String fullString = String.format("%.4f,%.4f", traitA, traitB);
		return fullString;
	}
	
	

}
