package ibid;

import java.util.ArrayList;

import cern.colt.matrix.DoubleMatrix2D;

public class AntigenicTree {
    
	//private static DoubleMatrix2D distMatrix;
    private static float[][] distMatrix;
	private static ArrayList<ArrayList<Float>> distArray;

	public AntigenicTree() {
		distMatrix = new float[1000][1000];
	}
	
	public void setdistMatrix(float a, int i, int j) {
	  distMatrix[i][j] = a;	
	}
	//
	public void updatedistMatrix(float da, int i, int j) {
		//source i infect target j
		//set da for [i,j]
		setdistMatrix(da, i, j);
		int numRows = 10; //should be replaced by virusid
		float colArray[] = new float[1000];
		
		//copy [,i] plus da to [,j]
		for(int row = 0; row < numRows; row++)
		{
		    colArray[row] = distMatrix[row][i];
		}
		for(int row = 0; row < numRows; row++)
		{
		    distMatrix[row][j] = colArray[row] + da;
		}
		//copy [,j] to [j,]
		for(int row = 0; row < numRows; row++)
		{
		    distMatrix[j][row] = distMatrix[row][j];
		}
		
		//set [j,j] = 0
		setdistMatrix(0, j ,j);
	}
}

