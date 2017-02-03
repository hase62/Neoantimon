package utility;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.analysis.SimpsonIntegrator;
import org.apache.commons.math.analysis.UnivariateRealFunction;

public class ExpRegression {
	private List <Double> X;
	   private List <Double> Y;
	   private double  powerOfX;
	   private List <Double> powX;
	   private List <Double> logY;
	   private LinearRegression L;
	   private int n;
	   private  double a;
	   private double b;
	   public ExpRegression(){
		   powerOfX = 1;
	   }
	   public void setPowerOfX(double d){
		   powerOfX = d;
		   if(!powX.isEmpty() && !X.isEmpty()){
			   powX.clear();
			   for(int i = 0; i < n ; i++){
				   powX.add(Math.pow(X.get(i), powerOfX));
			   }
			   
		   }
		   
	   }
	   public double getPowerOfX(){
		   return powerOfX;
	   }
	   public void setTrainingData(List <Double> X, List <Double> Y){
		   this.X = new ArrayList<Double>(X);
			this.Y = new ArrayList<Double>(Y);
		   n= Math.min(X.size(), Y.size());
		   powX = new ArrayList<Double>();
		   for(int i = 0; i < n ; i++){
			   powX.add(Math.pow(X.get(i), powerOfX));
		   }
		   logY = new ArrayList<Double>();
		   for(int i = 0; i < n ; i++){
			   logY.add(Math.log(Y.get(i)));
		   }
	   }
	   public void learn(List <Double> X, List <Double> Y){
		   setTrainingData(X,Y);
		   learn();
	   }
	   public void learn(){
		   L = new LinearRegression();
		   L.learn(powX, logY);
		   a = Math.exp(L.getIntercept());
		   b = L.getSlope();
	   }
	   public double predict(double x){
		   return a * Math.exp(b * Math.pow(x, powerOfX));
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
	   public void optimizePowerOfX(List <Double> powers) throws MyException{
		   Double minResidual = null;
		   LinearRegression minL = null;
		   Double bestPowerOfX = null;
		   double minPower = MyFunc.min(powers);
			double maxPower = MyFunc.max(powers);		   
		   int i;
		   for(i = 0; i < powers.size(); i++ ){
			   setPowerOfX(powers.get(i));
			   learn();
			   double r = getResidual();
			   if(minResidual == null || r < minResidual){
				   minResidual = r;
				   minL = L;
				   bestPowerOfX = powerOfX;
			   }	   
		   }
		   powerOfX = bestPowerOfX;
		   L = minL;
		   a = Math.exp(L.getIntercept());
		   b = L.getSlope();
		   if(getPowerOfX() == maxPower || getPowerOfX() == minPower){
				throw new MyException("PowerOfX (" + getPowerOfX() + ") might not be good!");
			}
	   }
	   
	   class Func implements UnivariateRealFunction{
		   public double value(double x) {
			   return predict(x);
		   }
	   }
	   public Double integrate(double LowerLimit, double upperLimit) throws MyException{
		   Func F = new Func();
		   SimpsonIntegrator SI = new SimpsonIntegrator(F);
		   try{
		   return SI.integrate(LowerLimit, upperLimit);
		   }
		   catch (Exception e) {
		     throw new MyException("integration has failed!");
		   }
	   }
	   public String toString(){
		   String S = "a = " + a + ", b = " + b + ", powerOfX = " + powerOfX + "\n";  
		   S += "X: " + MyFunc.join(" ", X) + "\n";
		   S += "Y: " + MyFunc.join(" ", Y) + "\n";
		   S += "powX: " + MyFunc.join(" ", powX) + "\n";
		   S += "logY: " + MyFunc.join(" ", logY) + "\n";
		   return S;
	   }
}
