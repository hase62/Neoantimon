package utility;

import java.util.*;
import org.apache.commons.cli.*;
import org.apache.commons.math.distribution.*;
import org.apache.commons.math.stat.inference.*;

import sun.reflect.Reflection;

public class MatrixCorrelation2{
	
	MyMat M1;
	MyMat M2;
 
	MyMat P;
	MyMat C;
 
	int minListSize = 6;
	double pcutoff = 10;
	double qcutoff = 10;
	boolean onlyMatchId = false;	
	boolean rawStatistics = false;
	
	
	public MatrixCorrelation2(MyMat m1, MyMat m2) throws Exception{		
		if(MyFunc.isect(m1.getColNames(), m2.getColNames()).size() >= minListSize){
			M1 = m1;
			M2 = m2;
		}else if(MyFunc.isect(m1.getRowNames(), m2.getColNames()).size() >= minListSize){
			M1 = m1;
			M2 = m2;
			M1.transpose();
		}else if(MyFunc.isect(m1.getColNames(), m2.getRowNames()).size() >= minListSize){
			M1 = m1;
			M2 = m2;
			M2.transpose();
		}else if(MyFunc.isect(m1.getRowNames(), m2.getRowNames()).size() >= minListSize){
			M1 = m1;
			M2 = m2;
			M1.transpose();
			M2.transpose();
		}else{
			throw new Exception("ERR: no row and column ids");
		}
		List <String> tmp = MyFunc.isect(M1.getColNames(), M2.getColNames());
		M1 = M1.getSubMatByCol(tmp);
		M2 = M2.getSubMatByCol(tmp);
		P = new MyMat(M1.getRowNames(), M2.getRowNames());
		C = new MyMat(M1.getRowNames(), M2.getRowNames());
	
	}
	
	
	public void calculatePvalue(){
		for(String s: M1.getRowNames()){
			for(String t: M2.getRowNames()){
				P.set(s, t, 1);
			}	
		}
		for(String s: M1.getRowNames()){
			L:for(String t: M2.getRowNames()){
				if(onlyMatchId & !s.equals(t)){
					continue L;
				}
				
				List <Double> L1 = M1.getRow(s);
				List <Double> L2 = M2.getRow(t);
				rmNan(L1, L2); 
				
				int n1 = (new HashSet<Double>(L1)).size();
				int n2 = (new HashSet<Double>(L2)).size();		
				
				double p = 1;
				if(n1 > 1 & n2 >1 & L1.size() > minListSize){
					try{
						if(n1==2 && n2==2){
							p = fisherExactTest(L1, L2);
						}else if(n1==2){
							p = tTest(L2, L1);		
						}else if(n2==2){
							p = tTest(L1, L2);
						}else {
							p = pearsonCorTest(L1, L2);
						}
					
						if(Double.isNaN(p)){
							System.err.println("WARN: unable to calculate a p value for " + s + " and " +t + " (NaN)");
						}
					
					}catch(Exception e){
						System.err.println("WARN: unable to calculate a p value for " + s + " and " +t);
					}
				}
				P.set(s, t, p);
				C.set(s, t, MyFunc.pearsonCorrelation(L1, L2)>0?1:-1);
			}
		}
	}
	
	
	
	private void rmNan (List <Double> L1, List <Double> L2){
		List <Double> l1 = new ArrayList<Double>(L1);
		List <Double> l2 = new ArrayList<Double>(L2);
		L1.clear();
		L2.clear();
		for(int i = 0; i < l1.size(); i++){
			if(!( Double.isNaN(l1.get(i)) ||   Double.isNaN(l2.get(i)))){
				L1.add(l1.get(i));
				L2.add(l2.get(i));	
			}
		}	
	}
	
	
	private double tTest (List <Double> continousL, List <Double> binaryL) throws Exception{
		List<Double> tmp = new ArrayList<Double>(new TreeSet<Double>(binaryL));
		List<Double> V = new ArrayList<Double>();
		List<Double> U = new ArrayList<Double>();
		for(int i = 0; i < continousL.size(); i++){
				if(eq(binaryL.get(i),tmp.get(0))){
					V.add(continousL.get(i));
				}else{
					U.add(continousL.get(i));
				}	
		}
		if(V.size()<=2 || U.size()<=2){
			throw new Exception();
		}
		
		double [] v = new double [V.size()];
		double [] u = new double [U.size()];
		for(int i = 0; i < v.length; i++){
			v[i] = V.get(i);
		}
		for(int i = 0; i < u.length; i++){
			u[i] = U.get(i);
		}
		
  	  TTestImpl TT = new  TTestImpl();
  	  return rawStatistics?Math.abs(TT.homoscedasticT(v,u)):TT.homoscedasticTTest(v,u);		
	}
	
	private boolean eq(double a, double b){
		return Math.abs(a-b)<0.0000000001;
	}
	
	private double pearsonCorTest(List<Double> L1, List<Double> L2) throws Exception{
			double r = MyFunc.pearsonCorrelation(L1, L2);
			int n = L1.size();
			if(eq(r,1)){
				return 0;
			}
			double t0 = Math.abs(r)*Math.sqrt(n-2)/Math.sqrt(1-Math.pow(r,2));
			if(rawStatistics){
				return Math.abs(t0);
			}else{
				TDistributionImpl T = new TDistributionImpl(n-2);
				
				return 2*(1-T.cumulativeProbability(t0));
			}
	}
	
   

	
	private double fisherExactTest(List<Double> L1, List<Double> L2) throws Exception{
		List<Double> tmp = new ArrayList<Double>(new TreeSet<Double>(L1));
		List<Double> tmp2 = new ArrayList<Double>(new TreeSet<Double>(L2));
		int a=0, b=0, c=0, d=0;
		for(int i = 0; i < L1.size(); i++){
			if(eq(L1.get(i), tmp.get(0)) & eq(L2.get(i), tmp2.get(0))){
				a++;
			}else if(eq(L1.get(i),tmp.get(0)) & eq(L2.get(i),tmp2.get(1))){
				b++;
			}else if(eq(L1.get(i), tmp.get(1)) & eq(L2.get(i), tmp2.get(0))){
				c++;
			}else if(eq(L1.get(i), tmp.get(1)) & eq(L2.get(i), tmp2.get(1))){
				d++;
			}
		}	
		return MyFunc.calculateFisheExactPvalue(a, b, c, d);
	}
	
	
	
	private void  printResults(){
		Map <String, Double> Pmap = new HashMap<String, Double>();
		Map <String, Double> Cmap = new HashMap<String, Double>();
		for(String s: M1.getRowNames()){
			for(String t: M2.getRowNames()){
				if(onlyMatchId){
					if(!s.equals(t)){
						continue;
					}else{
						Pmap.put(s, P.get(s,t));
						Cmap.put(s, C.get(s,t));
					}
				}else{
					Pmap.put(s + "\t" + t, P.get(s,t));	
					Cmap.put(s + "\t" + t, C.get(s,t));	
				}
			}
		}
		
		List <String> tmp = MyFunc.sortKeysByAscendingOrderOfValues(Pmap);
		
		if(tmp.size() <= 10){
			for(String s: tmp){
					System.out.println(s + "\t" + Cmap.get(s) + "\t" +  Pmap.get(s));
			}
		}else{
			Map <String, Double> Qmap = MyFunc.calculateQvalue(Pmap);
			for(String s: tmp){
				if(Qmap.get(s) <= qcutoff & Pmap.get(s) <= pcutoff ){
					System.out.println(s + "\t" + Cmap.get(s)+ "\t" + Pmap.get(s) + "\t" + Qmap.get(s));
				}
			}
		}
	}
 
	public MyMat getMinusLogPMatrix(){
		MyMat minusLogP = new MyMat(M1.getRowNames(), M2.getRowNames());
		if(rawStatistics){
			for(String s: M1.getRowNames()){
				for(String t: M2.getRowNames()){				
					minusLogP.set(s, t, C.get(s,t)*P.get(s,t));	
				}
			}
		}else{
			for(String s: M1.getRowNames()){
				for(String t: M2.getRowNames()){	
					double p;
					if(P.get(s,t)==0){
						p = Double.MAX_VALUE;
					}else if(P.get(s,t)==1){
						p = 0;
					}else{
						p = -Math.log10(P.get(s,t));
					}
					minusLogP.set(s, t, p);	
				}
			}
			double max = 20;
			for(String s: M1.getRowNames()){
				for(String t: M2.getRowNames()){				
					minusLogP.set(s, t, C.get(s,t)*((minusLogP.get(s,t) > max)?max:minusLogP.get(s,t)));	
				}
			}
		}
		return minusLogP;
	}
 
	private void printMinusLogPMatrix(){
	 System.out.print(getMinusLogPMatrix());
	}

 
	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("m", "pmat", false, "get minus log p-value matrix");
		options.addOption("p", "pcutoff", true, "cutoff for p-value");
		options.addOption("q", "pcutoff", true, "cutoff for q-value");
		options.addOption("C", "column", false, "get cor between columns (when the input is one matrix)");
		options.addOption("M", "match", false, "get cor only for matched IDs");
		options.addOption("s", "stat", false, "get raw statistics matrix");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabFle (tabFile)", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2 && argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabFile (tabFile)", options);
			return;
		}
			
		
		MatrixCorrelation2 MC;
		if(argList.size() == 2){
			MC = new MatrixCorrelation2(new MyMat(argList.get(0)), new MyMat(argList.get(1)));
		}else{
			MyMat M = new MyMat(argList.get(0));
			if(commandLine.hasOption("C")){
				M.transpose();
			}
			MC = new MatrixCorrelation2(M, M);
		}
		if(commandLine.hasOption("M")){
			MC.onlyMatchId = true;
		}
		if(commandLine.hasOption("s")){
			MC.rawStatistics = true;
		}
		if(commandLine.hasOption("p")){
			MC.pcutoff = Double.valueOf(commandLine.getOptionValue("p"));
		}
		if(commandLine.hasOption("q")){
			MC.qcutoff = Double.valueOf(commandLine.getOptionValue("q"));
		}
		MC.calculatePvalue();
		if(commandLine.hasOption("m") || commandLine.hasOption("s")){
			MC.printMinusLogPMatrix();
		}else{
			MC.printResults();
		}
	}
	
	
	
	
	
	
	

}
