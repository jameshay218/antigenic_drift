package ibid;
/* Holds random number genator necessities */
/* Trying to encapsulate this, so the RNG particulars can be changed if necessary */ 
/* Completely static class, allows no instances to be instantiated */

//import cern.jet.random.*;
import java.util.*;

public class Random {
								
	// methods

	public static int nextInt(int from, int to) {
		return cern.jet.random.Uniform.staticNextIntFromTo(from, to);
	}	
	
	public static double nextDouble() {
		return cern.jet.random.Uniform.staticNextDouble();		
	}
	
	public static double nextDouble(double from, double to) {
		return cern.jet.random.Uniform.staticNextDoubleFromTo(from, to);		
	}	

	public static double nextNormal() {
		return cern.jet.random.Normal.staticNextDouble(0.0,1.0);
	}
	
	public static double nextNormal(double mean, double sd) {
		return cern.jet.random.Normal.staticNextDouble(mean,sd);
	}	

	// tuned with mean
	public static double nextExponential(double lambda) {
		return cern.jet.random.Exponential.staticNextDouble(1.0/lambda);
	}
	
	// tuned with alpha and beta, matching Mathematica's notation
	public static double nextGamma(double alpha, double beta) {
		return cern.jet.random.Gamma.staticNextDouble(alpha, 1/beta);
	}
        
        // tuned with alpha and beta, matching Mathematica's notation
	public static double nextColtGamma(double alpha, double beta) {
		return cern.jet.random.Gamma.staticNextDouble(alpha, beta);
	}
	
	public static int nextPoisson(double lambda) {
		return cern.jet.random.Poisson.staticNextInt(lambda);
	}
	
	public static boolean nextBoolean(double p) {
		boolean x = false;
		if (nextDouble() < p) {
			x = true;
		}
		return x;
	}	
	

}
