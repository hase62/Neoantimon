package tensor;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.Options;
import org.apache.commons.math.random.JDKRandomGenerator;
import org.apache.commons.math.random.RandomGenerator;

import utility.*;

public class Order3TensorGenerator {
	int n1 = 300;
	int n2 = 1000;
	int n3 = 20;
	
	double subTensorRatio = 0.3;
	double subTensorNumber = 10;
	double signalStrength = 3;
	
	
	RandomGenerator RG;
	Order3Tensor T;
	
	public Order3TensorGenerator() {
		RG = new JDKRandomGenerator();
		T = new Order3Tensor(n1, n2, n3);
	}
	
	protected void simulate(){
		int i,j,k,l;
		for(i = 0; i < n1; i++){
			for(j = 0; j < n2; j++){
				for(k = 0; k < n3; k++){
					T.set(i,j,k, RG.nextGaussian());
				}
			}
		}
		
		for(l=0;l<subTensorNumber;l++){
			List <String> order1cluster  = new ArrayList<String>(MyFunc.sample(T.getName1(), (int)Math.round(subTensorRatio*n1)));
			List <String> order2cluster  = new ArrayList<String>(MyFunc.sample(T.getName2(), (int)Math.round(subTensorRatio*n2)));
			List <String> order3cluster  = new ArrayList<String>(MyFunc.sample(T.getName3(), (int)Math.round(subTensorRatio*n3)));
			
			double sign = RG.nextGaussian()>0?1:-1;
			for(String s1: order1cluster){
				for(String s2: order2cluster){
					for(String s3: order3cluster){
						T.set(s1,s2,s3, sign*signalStrength + RG.nextGaussian());
					}
				}
			}
		}
	}

	
	public static void main(String [] args) throws Exception{
		Order3TensorGenerator G = new Order3TensorGenerator();
		G.simulate();
		//System.out.print(G.T);
	
		ClusteredOrder3Tensor T  = new ClusteredOrder3Tensor(G.T);
		T.performClustering();
		System.out.print(T);
	}
	
	

}
