package utility;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.DataFormatException;

import sun.reflect.Reflection;

import weka.classifiers.functions.SimpleLinearRegression;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;

public class LinearRegression {
	
	 private SimpleLinearRegression L;
	   private Instances insts;
	   private Attribute atX;
	   private Attribute atY;
	   private Double slope;
	   private Double intercept;
	   private List <Double>  X;
	   private List <Double> Y;
	   private int n;
	   
	   public LinearRegression() {
		   atX = new Attribute("X");
			atY = new Attribute("Y");
			FastVector att = new FastVector();
			att.addElement(atX);
			att.addElement(atY);
			insts = new Instances("tmp",att,3);
			insts.setClass(atY);
			L = new SimpleLinearRegression();
			slope = null;
			intercept = null;
			X = null;
			Y = null;
	   }
	   public LinearRegression(List <Double> X, List<Double> Y){
		   atX = new Attribute("X");
			atY = new Attribute("Y");
			FastVector att = new FastVector();
			att.addElement(atX);
			att.addElement(atY);
			insts = new Instances("tmp",att,3);
			insts.setClass(atY);
			L = new SimpleLinearRegression();
			slope = null;
			intercept = null;
			this.X = null;
			this.Y = null;
			setTrainingData(X,Y);
	   }
	   public void setTrainingData(List <Double> X, List<Double> Y){
		   this.X = new ArrayList<Double>(X);
			this.Y = new ArrayList<Double>(Y);
		   n= Math.min(X.size(), Y.size());
	   }

	   public void learn(List <Double> X, List<Double> Y) {
		   setTrainingData(X,Y);
		   learn();
	   }
		public void learn() {   
		   int i;
		   for(i=0; i<n ;i++){
			   Instance inst = new Instance(2);
			   inst.setValue(atX,X.get(i));
			   inst.setValue(atY,Y.get(i));
			   insts.add(inst);
		   }
		   try{
		   L.buildClassifier(insts);
		   }
		   catch (Exception e) {
			}
		   slope = L.getSlope();
		   intercept = L.getIntercept();
	   }
	   public Double getSlope(){
		   return slope;
	   }
	   public Double getIntercept(){
		   return intercept;
		   
	   }
	   public Double predict(double x){
		   return slope * x + intercept;
	   }
	   public List <Double> predict(List<Double> X){
		   List <Double> Y = new ArrayList<Double>();
		   int i;
		   for(i = 0; i< X.size(); i++){
			   Y.add(predict(X.get(i)));
		   	}
		   return Y;   
	   }
	   public double getResidual(){
		   List <Double> preY = predict(X);
		   List <Double> tmp = new ArrayList<Double>();
		   int i;
		   for(i = 0; i < n; i++){	 
			   tmp.add(Math.pow((preY.get(i)- Y.get(i)),2));
		   }
		  return MyFunc.sum(tmp); 
	   }
	   public String toString(){
		   String S = "slope = " + slope + ", intercept = " + intercept;  
		   return S;
	   }
	   
	   public static void main(String[] args) throws Exception {
		   if(args.length != 2){
				System.err.println(Reflection.getCallerClass( 1 ).getName() + " training_data test_data");
				return;
			}
		   BufferedReader inputStream = new BufferedReader(new FileReader(args[0]));
			String line;
	 		List<Double> trainingX = new ArrayList<Double>();
	 		List<Double> trainingY = new ArrayList<Double>();
	 		while((line = inputStream.readLine()) != null){
	 			List<String> str = Arrays.asList(line.split("\t"));
				if(str.size() < 2){
					throw new DataFormatException("LinearRegression: format is wrong!");
				}
				trainingX.add(Double.valueOf(str.get(0)));
				trainingY.add(Double.valueOf(str.get(1)));
			}
	 		LinearRegression L = new LinearRegression(trainingX, trainingY);
	 		L.learn();
	 		inputStream = new BufferedReader(new FileReader(args[1]));
	 		while((line = inputStream.readLine()) != null){
	 			List<String> str = Arrays.asList(line.split("\t"));
				if(str.size() == 0){
					throw new DataFormatException("LinearRegression: format is wrong!");
				}
				double x = Double.valueOf(str.get(0));
				System.out.println(x + "\t" + L.predict(x));
			}
	   }
	
	
	
	
	
	
	

}
