package ibid;

//import Phenotype;
//import Virus;

public class Virus {

	// fields
	private Virus parent;		// parent virus 
	private double birth;		// measured in years relative to burnin
	private int cluster;		// antigenic type
	private int deme;
	
	public Phenotype phenotype;

	public Virus() {
		phenotype = new Phenotype();
	}
	
	public Virus(Virus p) {
		setParent(p);
		this.phenotype = p.getPhenotype();
	}
	
    private void setParent(Virus p) {
    	parent = p;
    }
    
    public Phenotype getPhenotype() {
    	return new Phenotype(phenotype);
    }
    
    public void mutate() {

    }
}
