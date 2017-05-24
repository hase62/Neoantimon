package bct;

import java.util.*;

import org.apache.commons.cli.*;

import sun.reflect.Reflection;
import utility.*;

public class BCTcorrelationTest {
	BCT M1;
	BCT M2;
	
	List <String> r1;
	List <String> r2;
	List <String> c;
	
	
	MyMat P1;
	MyMat P2;
	
	int n1;
	int n2;
	int m;

	double posWeight = 1;
	double sampWeight = 1;
	
	boolean comb = false;
	boolean larger = true;
	boolean intercor = false;
	
	double replaceNan = -1;
	
	Map<String, Double> Pvalue;
	
	PoibinRNA PB;
	
	double pcutoff = 3;
	
	public static class Operation{
		private final String name;
		private Operation(String name){this.name = name;}
		public Operation(Operation b){this.name = b.name;}
		@Override
		public String toString(){return name;}
		public boolean equals(Operation b){return name.equals(b.toString());};
		public static final Operation  AND = new Operation("AND");
		public static final Operation OR = new Operation("OR");
		public static final Operation XOR = new Operation("XOR");
	}
	
	
	interface PBstatistic{
		int get(int i, int j);
	}
	
	interface PBparameter{
		List<Double> get(int i, int  j);
	}
	
	interface Pcalculation{
		double get(int i, int  j);
	}
	
	
	
	Operation operation;
	PBstatistic statistic;
	PBparameter parameter;
	Pcalculation pcalculation;
	
	private void setAND(){
		operation =  Operation.AND;
		statistic = new PBstatistic(){
			public int get(int i, int j){
				int d = 0;
				for(int k  = 0; k < m;k++){
					if(M1.is1(i,k) & M2.is1(j,k)){
						d++;
					}
				}
				return d;
			}
		};
		
		parameter = new PBparameter(){
			public List<Double> get(int i, int j){
				List <Double> p = new ArrayList<Double>();
				for(int k  = 0; k < m;k++){
					p.add(P1.get(i,k)*P2.get(j,k));
				}
				return p;
			}
		};
	}
	
	private void setOR(){
		operation =  Operation.OR;
		statistic = new PBstatistic(){
			public int get(int i, int j){
				int d = 0;
				for(int k  = 0; k < m;k++){
					if(M1.is1(i,k) | M2.is1(j,k)){
						d++;
					}
				}
				return d;
			}
		};
		
		parameter = new PBparameter(){
			public List<Double> get(int i, int j){
				List <Double> p = new ArrayList<Double>();
				for(int k  = 0; k < m;k++){
					p.add(1- (1-P1.get(i,k))*(1 -P2.get(j,k)));
				}
				return p;
			}
		};
		
	}
	
	private void setXOR(){
		operation =  Operation.XOR;
		statistic = new PBstatistic(){
			public int get(int i, int j){
				int d = 0;
				for(int k  = 0; k < m;k++){
					if((M1.is1(i,k) & M2.is0(j,k) ) | (M1.is0(i,k) & M2.is1(j,k))){
						d++;
					}
				}
				return d;
			}
		};
		
		parameter = new PBparameter(){
			public List<Double> get(int i, int j){
				List <Double> p = new ArrayList<Double>();
				for(int k  = 0; k < m;k++){
					p.add(P1.get(i,k)*(1-P2.get(j,k)) + (1-P1.get(i,k))*P2.get(j,k));
				}
				return p;
			}
		};
	}
	
	private void testLarger(){
		larger = true;
		pcalculation  = new Pcalculation(){
			public double  get(int i, int j){
				int k =  statistic.get(i, j);
				List <Double> p = parameter.get(i,j);
				return 1- PB.getCdf(k,p);
			}
		};
	}
	
	private void testSmaller(){
		larger = false;
		pcalculation  = new Pcalculation(){
			public double  get(int i, int j){
				int k =  statistic.get(i, j);
				List <Double> p = parameter.get(i,j);
				return PB.getCdf(k-1,p);
			}
		};
	}
	
	private void testPcutoff(double d){
		pcutoff = d;
	}
	
	public  void setPosWeight(double d){
		if(d <= 1 & d >= 0){
			posWeight = d;
		}
	}
	
	public  void setSampWeight(double d){
		if(d <= 1 & d >= 0){
			sampWeight = d;
		}
	}
	
	public BCTcorrelationTest(MyMat A, MyMat B){
		c = MyFunc.isect(A.getColNames(), B.getColNames());
		r1 = A.getRowNames();
		r2 = A.getRowNames();
		M1 = new BCT(A.getSubMatrix(r1, c));
		M2= new BCT(B.getSubMatrix(r2, c));
		n1 = M1.rowSize();
		n2 = M2.rowSize();
		m = c.size();
		PB = new PoibinRNA();
		setAND();
		testLarger();
		comb=true;
		intercor=false;
	}
	
	public BCTcorrelationTest(MyMat A){
		M1 = new BCT(A);
		M2 = new BCT(A);
		n1 = M1.rowSize();
		n2 = n1;
		m = M1.colSize();
		PB = new PoibinRNA();
		setAND();
		testLarger();
		comb=true;
		intercor=true;
	}
	
	public void testMatchedRows(){
		if(!intercor){
			comb = false;
			List<String> r = MyFunc.isect(r1, r2);
			
			List <Integer> R1 = new ArrayList <Integer>();
			for(int i = 0; i < r1.size(); i++){
				if(r.contains(r1.get(i))){
					R1.add(i);
				}
			}
			
			List <Integer> R2 = new ArrayList <Integer>();
			for(int i = 0; i < r2.size(); i++){
				if(r.contains(r2.get(i))){
					R2.add(i);
				}
			}

			r1=r;
			r2=r;
			
			List <Integer> C = new ArrayList <Integer>();
			for(int i = 0; i < c.size(); i++){
				C.add(i);
			}
			
			M1 = M1.getSubTable(R1, C);
			M2= M2.getSubTable(R2, C);			
			n1 = M1.rowSize();
			n2 = M2.rowSize();
		}
	}
	
	
	private static List <Double> reWeight (List <Double> L, double k){
		List <Double> tmp = new ArrayList <Double>();
		List <Double> tmp2 = new ArrayList <Double>();
		for(double d: L){
			tmp.add(Math.pow(d, k));
		}
		double tmp3 = MyFunc.sum(tmp);
		for(double d: tmp){
			tmp2.add(d/tmp3);
		}
		return tmp2;
	}
	
	private MyMat calculateParameters(BCT M){
		double a;
		List <Double> q;
		List <Double> r;
		a = M.sum();
		q = new ArrayList<Double>();
		r = new ArrayList<Double>();
		
		for(int i = 0; i < M.rowSize(); i++){
			q.add(M.rowSum(i)/a);
		}
		
		if(posWeight  != 1){
			q = reWeight(q,posWeight);
		}	
		
		
		for(int k = 0; k < M.colSize(); k++){
			r.add(M.colSum(k)/a);
		}
		if(sampWeight != 1){
			r = reWeight(r,sampWeight);
		}
		MyMat P = new MyMat(M.rowSize(), M.colSize());
		for(int i = 0; i < M.rowSize(); i++){
			for(int k = 0;k < M.colSize(); k++){
				double p = a*q.get(i)*r.get(k);
				if(p > 1){
					System.err.println(i + " " + k + " "+ p);
				}
				P.set(i, k, p);
			}
		}
		return P;
	}
	
	
	public void calculateParameters(){
		P1 = calculateParameters(M1);
		P2 = calculateParameters(M2);
	}
	
	public void getPvalue(){
		if(!comb){
			Pvalue = new HashMap<String, Double>();
			for(int i = 0; i<n1; i++){
				Pvalue.put(r1.get(i), pcalculation.get(i,i));
			}		
		}else{
			Pvalue = new HashMap<String, Double>();
			for(int i = 0; i<n1; i++){
				for(int j = 0; j<n2; j++){
					if(intercor & (j >= i)){
						continue;
					}
					Pvalue.put(r1.get(i)+ " "+r2.get(j), pcalculation.get(i,j));
				}		
			}
		}
	}
	
	public void convertPvalue2LogScale(){
		double  tmp = 1;
		for(String p: Pvalue.keySet()){
			if(Pvalue.get(p) != 0 &  Pvalue.get(p) < tmp){
				tmp = Pvalue.get(p);
			}
		}
		double max = - Math.log10(tmp);
		
		
		for(String p: Pvalue.keySet()){
			tmp = Pvalue.get(p);
			if(tmp==1){
				tmp=0;
			}else if(tmp==0){
				tmp = max;
			}else if(Double.isNaN(tmp)){
				tmp = replaceNan;
			}else{
				tmp = - Math.log10(tmp);
			}
			
			
			Pvalue.put(p, tmp);
		}
		
		Map <String, Double> tmp2 = new HashMap<String, Double>(Pvalue);
		Pvalue = new LinkedHashMap<String, Double>();
		for(String s: MyFunc.sortKeysByDescendingOrderOfValues(tmp2)){
			if(tmp2.get(s) > pcutoff){
				Pvalue.put(s, tmp2.get(s));
			}
		}
	}
	 
	
	public void perform(){
		calculateParameters();
		getPvalue();
		convertPvalue2LogScale();
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("o", "or", false, "test OR");
		options.addOption("x", "eor", false, "test XOR");
		options.addOption("s", "small", false, "test smaller");
		options.addOption("p", "pcutoff", true, "pvalue cutoff");
		options.addOption("m", "match", false, "test matched rows");
		options.addOption("P", "suppos", true, "suppress position-wise weight");
		options.addOption("S", "samppos", true, "suppress sample-wise weight");
		
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] input.tab", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(!(argList.size() == 1 | argList.size() == 2)){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] input.tab", options);
			System.exit(1);
		}
		
		BCTcorrelationTest B = null;
		if(argList.size() == 1 ){
			B = new BCTcorrelationTest(new MyMat(argList.get(0)));
		}else{
			B = new BCTcorrelationTest(new MyMat(argList.get(0)), new MyMat(argList.get(1)));
		}
	
		if(commandLine.hasOption("o")){
			B.setOR();
		}
		if(commandLine.hasOption("x")){
			B.setXOR();
		}
		if(commandLine.hasOption("s")){
			B.testSmaller();
		}
		if(commandLine.hasOption("p")){
			B.testPcutoff(Double.valueOf(commandLine.getOptionValue("p")));
		}
		if(commandLine.hasOption("m")){
			B.testMatchedRows();
		}
		if(commandLine.hasOption("S")){
			B.setSampWeight(Double.valueOf(commandLine.getOptionValue("S")));
		}
		if(commandLine.hasOption("P")){
			B.setPosWeight(Double.valueOf(commandLine.getOptionValue("P")));
		}
		
		B.perform();
		
		for(String s: B.Pvalue.keySet()){
			System.out.println(s + "\t" + B.Pvalue.get(s));
		}
	}

}
