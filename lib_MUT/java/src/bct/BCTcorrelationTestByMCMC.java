package bct;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;
import utility.MyFunc;
import utility.MyMat;

public class BCTcorrelationTestByMCMC {
	BCT M1;
	BCT M2;
	
	BCTsampler Bsamp;
	BCTsampler Bsamp2;
	int NforBurnIn = 1000;
	int N = 100000;
	
	List <String> r1;
	List <String> r2;
	List <String> c;
		
	int n1;
	int n2;
	int m;
	
	boolean comb = false;
	boolean larger = true;
	boolean intercor = false;
	
	double replaceNan = -1;
	
	Map<String, Double> Statistic;
	Map<String, Double> Pvalue;
	
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
		double get(int i, int j);
		double getFromBsamp(int i, int j);
	}
	
	Operation operation;
	PBstatistic statistic;
	
	private void setAND(){
		operation =  Operation.AND;
		statistic = new PBstatistic(){
			public double get(int i, int j){
				double d = 0;
				for(int k  = 0; k < m;k++){
					if(M1.is1(i,k) & M2.is1(j,k)){
						d++;
					}
				}
				return d;
			}
			
			public double getFromBsamp(int i, int j){
				double d = 0;
				for(int k  = 0; k < m;k++){
					if(Bsamp.getCurrentBCT().is1(i,k) & Bsamp2.getCurrentBCT().is1(j,k)){
						d++;
					}
				}
				return d;
			}
			
		};
	}
		
	private void setOR(){
		operation =  Operation.OR;
		statistic = new PBstatistic(){
			public double get(int i, int j){
				double d = 0;
				for(int k  = 0; k < m;k++){
					if(M1.is1(i,k) | M2.is1(j,k)){
						d++;
					}
				}
				return d;
			}
			public double getFromBsamp(int i, int j){
				double d = 0;
				for(int k  = 0; k < m;k++){
					if(Bsamp.getCurrentBCT().is1(i,k) | Bsamp2.getCurrentBCT().is1(j,k)){
						d++;
					}
				}
				return d;
			}
		};		
	}
		
	private void setXOR(){
		operation =  Operation.XOR;
		statistic = new PBstatistic(){
			public double get(int i, int j){
				double d = 0;
				for(int k  = 0; k < m;k++){
					if((M1.is1(i,k) & M2.is0(j,k) ) | (M1.is0(i,k) & M2.is1(j,k))){
						d++;
					}
				}
				return d;
			}
			public double getFromBsamp(int i, int j){
				double d = 0;
				for(int k  = 0; k < m;k++){
					if((Bsamp.getCurrentBCT().is1(i,k) & Bsamp2.getCurrentBCT().is0(j,k) ) | (Bsamp.getCurrentBCT().is0(i,k) & Bsamp2.getCurrentBCT().is1(j,k))){
						d++;
					}
				}
				return d;
			}
		};	
	}
		
	private void testLarger(){
		larger = true;
	}
		
	private void testSmaller(){
		larger = false;
	}
		
	private void setPcutoff(double d){
		pcutoff = d;
	}
	
	private void setN(int i){
		N = i;
	}
	private void setNforBurnIn(int i){
		NforBurnIn = i;
	}
	
		
		
	public BCTcorrelationTestByMCMC(MyMat A, MyMat B){
		c = MyFunc.isect(A.getColNames(), B.getColNames());
		r1 = A.getRowNames();
		r2 = B.getRowNames();
		M1 = new BCT(A.getSubMatrix(r1, c));
		M2= new BCT(B.getSubMatrix(r2, c));
		n1 = M1.rowSize();
		n2 = M2.rowSize();
		m = c.size();
		setAND();
		testLarger();
		comb=true;
		intercor=false;
		Bsamp = new BCTsampler(M1);
		Bsamp2 = new BCTsampler(M2);
	}
		
	public BCTcorrelationTestByMCMC(MyMat A){
		c = A.getColNames();
		r1 = A.getRowNames();
		r2 = A.getRowNames();
		M1 = new BCT(A);
		M2 = M1;
		n1 = M1.rowSize();
		n2 = n1;
		m = M1.colSize();
		setAND();
		testLarger();
		comb=true;
		intercor=true;
		Bsamp = new BCTsampler(M1);
		Bsamp2 = Bsamp;
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
		
	public void getStatistics(){
		Statistic = new HashMap<String, Double>();
		if(!comb){
			for(int i = 0; i<n1; i++){
				Statistic.put(r1.get(i), statistic.get(i, i));
			}
		}else{
			for(int i = 0; i<n1; i++){
				for(int j = 0; j<n2; j++){
					if(intercor & (j >= i)){
						continue;
					}
					Statistic.put(r1.get(i)+ " "+r2.get(j), statistic.get(i,j));
				}		
			}	
		}
	}
		
	public void getPvalue(){
		Pvalue = new HashMap<String, Double>();
		if(!comb){
			for(int i = 0; i<n1; i++){
				Pvalue.put(r1.get(i), 0.0);
			}		
			for(int c =0; c < NforBurnIn; c++){
				Bsamp.mcmc.next();
			}
		
			for(int c = 0; c < N; c++){
				Bsamp.mcmc.next();
				for(int i = 0; i<n1; i++){
					if(larger){
						if(statistic.getFromBsamp(i, i) > Statistic.get(r1.get(i))){
							Pvalue.put(r1.get(i), Pvalue.get(r1.get(i)+1));
						}
					}else{
						if(statistic.getFromBsamp(i, i) < Statistic.get(r1.get(i))){
							Pvalue.put(r1.get(i), Pvalue.get(r1.get(i)+1));
						}	
					}
				}
			}
			for(int i = 0; i<n1; i++){
				if(Pvalue.get(r1.get(i))==0.0){
					Pvalue.put(r1.get(i), 1.0);
				}
				Pvalue.put(r1.get(i), Pvalue.get(r1.get(i))/N);
			}	
		}else{
			for(int i = 0; i<n1; i++){
				for(int j = 0; j<n2; j++){
					if(intercor & (j >= i)){
						continue;
					}
					Pvalue.put(r1.get(i)+ " "+r2.get(j), 0.0);
				}		
			}
			for(int c =0; c < NforBurnIn; c++){
				Bsamp.mcmc.next();
				if(!intercor){
					Bsamp2.mcmc.next();
				}
			}
			for(int c = 0; c < N; c++){
				Bsamp.mcmc.next();
				if(!intercor){
					Bsamp2.mcmc.next();
				}
				for(int i = 0; i<n1; i++){
					for(int j = 0; j<n2; j++){
						if(intercor & (j >= i)){
							continue;
						}
						if(larger){
							if(statistic.getFromBsamp(i, j) > Statistic.get(r1.get(i)+ " "+r2.get(j))){
								Pvalue.put(r1.get(i)+ " "+r2.get(j), Pvalue.get(r1.get(i)+ " "+r2.get(j))+1);
							}
						}else{
							if(statistic.getFromBsamp(i, j) < Statistic.get(r1.get(i)+ " "+r2.get(j))){
								Pvalue.put(r1.get(i)+ " "+r2.get(j), Pvalue.get(r1.get(i)+ " "+r2.get(j))+1);
							}	
						}
					}
				}
			}
			
			for(int i = 0; i<n1; i++){
				for(int j = 0; j<n2; j++){
					if(intercor & (j >= i)){
						continue;
					}
					if(Pvalue.get(r1.get(i)+ " "+r2.get(j))==0.0){
						Pvalue.put(r1.get(i)+ " "+r2.get(j), 1.0);
					}
					Pvalue.put(r1.get(i)+ " "+r2.get(j), Pvalue.get(r1.get(i)+ " "+r2.get(j))/N);
					System.err.println(Pvalue.get(r1.get(i)+ " "+r2.get(j)));
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
		getStatistics();
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
		options.addOption("n", "nchain", true, "the number of chains");
		options.addOption("b", "nbin", true, "the number of burn in");
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
			
		BCTcorrelationTestByMCMC B = null;
		if(argList.size() == 1 ){
			B = new BCTcorrelationTestByMCMC(new MyMat(argList.get(0)));
		}else{
			B = new BCTcorrelationTestByMCMC(new MyMat(argList.get(0)), new MyMat(argList.get(1)));
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
			B.setPcutoff(Double.valueOf(commandLine.getOptionValue("p")));
		}
		if(commandLine.hasOption("n")){
			B.setN(Integer.valueOf(commandLine.getOptionValue("n")));
		}
		if(commandLine.hasOption("b")){
			B.setNforBurnIn(Integer.valueOf(commandLine.getOptionValue("b")));
		}
		if(commandLine.hasOption("m")){
			B.testMatchedRows();
		}			
		B.perform();
		for(String s: B.Pvalue.keySet()){
			System.out.println(s + "\t" + B.Pvalue.get(s) + "\t" + B.Statistic.get(s));
		}
	}
}

	

