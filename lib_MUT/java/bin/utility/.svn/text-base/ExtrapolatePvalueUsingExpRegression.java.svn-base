package utility;

import java.util.ArrayList;
import java.util.List;

public class ExtrapolatePvalueUsingExpRegression {
	
	
	List <Double> nullDist;
	double minPercentile = 0.7;	
	double maxPercentile = 0.9;
	double minPower = 1.0;
	double maxPower = 3.0;
	
	MyFunc.Density D;
	ExpRegression ER;
	
	double upperLimit;
	
	public  ExtrapolatePvalueUsingExpRegression (List <Double> nullDist){
		this.nullDist = nullDist;
		upperLimit = MyFunc.mean(nullDist)*50;
	}
	
	
	public void setXrange(double min, double max){
		 minPercentile = min;
		 maxPercentile = max;
	}
	
	public void setPoweRange(double min, double max){
		 minPower = min;
		 maxPower = max;
	}
	
	public void setUpperLimit(double d){
		upperLimit = d;
	}
	
	
	public void regress(){
		D = new MyFunc.Density(nullDist);
		double minX = MyFunc.percentile(nullDist, minPercentile);
		double maxX = MyFunc.percentile(nullDist, maxPercentile);
		double d = (maxX - minX) / 500;
		List <Double> X = new ArrayList<Double>();
		double x;
		for(x = minX; x <= maxX; x += d ){
			X.add(x);
		}
		List <Double> Y = D.estimate(X);
		ER  = new ExpRegression();
		ER.setTrainingData(X, Y);
		List <Double> powers = new ArrayList<Double>();
	
		for( double p = minPower; p <= maxPower; p += 0.1){
			powers.add(p);
		}
		try {
			ER.optimizePowerOfX(powers);
		} catch (MyException e) {
			e.printStackTrace();
		}
	}
	
	
	public double  estimatePvalue(double d){
		return  ER.integrate(d, upperLimit);
	}
	
	












}




