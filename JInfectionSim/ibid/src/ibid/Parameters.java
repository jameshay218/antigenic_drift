package ibid;

public class Parameters {

		// global parameters
		public static int day = 365*10;
		public static int N = 100000;		//total population size
		public static int initialS = 70000;
		public static int initialI = 700;
		public static int initialR = 29300;
		public static double mu = 1/(70*365); 	//birth and death rate
		public static double rho = 0.4;		//probability of successful transmission
		public static double contact = 1; 		//contact rate
		public static double wan = 1/(0.5*365);	 	//immune waning rate
		public static double gamma = 1/3.8;	//recovery rate

		public static double getBirth() {
			return mu;
		}
}
