package utility;

import java.io.*;
import java.util.*;

import org.apache.commons.math.random.*;
import org.apache.commons.math3.distribution.NormalDistribution;

public class PvalueCalculatorWithGPD {
	
	private double a;
	private double k;
	
	private double a_se;
	private double k_se;
	
	private int Nexc = 250;
	private int M = 10;
	private int N;
	private double t;
	
	boolean testLarger = true;
	
	List <Double> empDist;
	public PvalueCalculatorWithGPD(){}	
	public PvalueCalculatorWithGPD(List <Double> D){
		 setEmpiricalDistribution(D);
		 setT();
		 setGPDparameter();
	}
	
	public void setEmpiricalDistribution(List <Double> D){
		empDist = D;
		N = empDist.size();	
	}
	
	public void testLess(){
		testLarger = false;
	}
	
	public void setT(){
		Collections.sort(empDist); //ascending order
		Collections.reverse(empDist);
		t = (empDist.get(Nexc-1)+empDist.get(Nexc))*0.5;
	}
	
	public String toString(){
		return "a=" + a + "+-" + a_se  + ", k=" + k + "+-" + k_se + ", t=" + t; 
	}
	
	public void  setGPDparameter(){
		 String tmpFile = "tmp" + Math.round(Math.random()*100000000); 
		   try {
			   PrintWriter os;
			   os = new PrintWriter(new FileWriter(tmpFile + ".nulldist"));
			   for(int i = 0, n = empDist.size(); i < n; i++){
				   os.println(empDist.get(i).toString());
			   }
			   os.flush();
			   os.close();
			   os = new PrintWriter(new FileWriter(tmpFile + ".R"));
			   os.println(" data<-read.table(\"" + tmpFile + ".nulldist" + "\")[,1]");
			   os.println("library(\"ismev\")");
			   os.println("fit<-gpd.fit(data," + t + ")");
			   os.println("if(fit$conv==0){");
			   os.println("write(c(fit$mle, fit$se), \""+ tmpFile + ".mle" + "\", ncolumns=1)");
			   os.println("}");
			   os.flush();
			   os.close();
			   MyFunc.runRscript(tmpFile + ".R");
			   BufferedReader is = new BufferedReader(new FileReader( tmpFile + ".mle"));
			   a = Double.valueOf(is.readLine());
			   k = Double.valueOf(is.readLine());
			   a_se = Double.valueOf(is.readLine());
			   k_se = Double.valueOf(is.readLine());
			   is.close();
			   File f = new File(tmpFile + ".nulldist");
			   if(!f.delete()){
				   f.deleteOnExit();
			   }
			   f = new File(tmpFile + ".R");
			   if(!f.delete()){
				   f.deleteOnExit();
			   }
			   f = new File(tmpFile + ".mle");
			   if(!f.delete()){
				   f.deleteOnExit();
			   }
		   } catch (Exception e) {
			   //e.printStackTrace();
			   File f = new File(tmpFile + ".nulldist");
			   if(!f.delete()){
				   f.deleteOnExit();
			   }
			   f = new File(tmpFile + ".R");
			   if(!f.delete()){
				   f.deleteOnExit();
			   }
			   f = new File(tmpFile + ".mle");
			   if(!f.delete()){
				   f.deleteOnExit();
			   }
			   throw  new ArithmeticException();
		   }
	}
	
	private double CDF(double x){
		 double z = x-t;
		 if(k==0){
			 return 1- Math.exp(-z/a);
		 }else{
			 return 1- Math.pow((1+k*z/a),-1/k);
		 }
	}
	
	public double getPgpd(double x){
		return ((double)Nexc)/N*(1-CDF(x));
	}
	
	public double getPecdf(double x){
		int j = 0;
		for(int i = 0; i<N; i++){
			if(empDist.get(i) >= x){
				j++;
			}else{
				break;
			}
		}
		return (double)(1+j)/N;
	}
	public double getPvalue(double x){
		int j = 0;
		for(int i = 0; i<N; i++){
			if(empDist.get(i) >= x){
				j++;
			}else{
				break;
			}
		}
		
		if(j>=M){
			return getPecdf(x);
		}else{
			return getPgpd(x);
		}	
	}

	public static void main(String[] args) throws Exception {
		RandomGenerator rg = new JDKRandomGenerator();
		GaussianRandomGenerator GRG = new GaussianRandomGenerator(rg);
		List <Double> D = new ArrayList<Double>();
		for(int i = 0; i < 10000; i++){
			D.add(GRG.nextNormalizedDouble());
		}
		NormalDistribution ND = new NormalDistribution();
		PvalueCalculatorWithGPD P = new PvalueCalculatorWithGPD(D);
		System.err.println(P);
		for(double x = 0; x < 5; x+=0.5){
			double p = 1.0 - ND.cumulativeProbability(x);
			double p2 = P.getPvalue(x);
			double p3 = P.getPecdf(x);
			System.out.println(x + "\t" + p + "\t" + p2  + "\t" + p3);
		}
	
	}
	
}
