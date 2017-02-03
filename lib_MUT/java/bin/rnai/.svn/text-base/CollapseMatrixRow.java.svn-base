package rnai;

import sun.reflect.Reflection;
import utility.*;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;

import org.apache.commons.cli.*;

public class CollapseMatrixRow {
	MyMat R; //original
	MyMat S; // rank
	MyMat T; // statistics
	MyMat P; // p-value or z-score
	Map <String, List<String>> gene2probe;
	int nullDistSize = 100000;
	boolean ascending = true;
	boolean zscore = false;
	int  rowSize;
	Map <Integer, List<Double> > nullDist;
	List <Integer> sampledIndex;
	
	public  CollapseMatrixRow(MyMat r, Map <String, String> probe2gene) throws Exception{
		R = r;
		List <String> tmp = new ArrayList<String>(probe2gene.keySet());
		List <String> tmp2 = MyFunc.isect(tmp, r.getRowNames());
		R = R.getSubMatByRow(tmp2);
		gene2probe = new HashMap <String, List<String>> ();
		for(String s: tmp2){
			if(!gene2probe.containsKey(probe2gene.get(s))){
				gene2probe.put(probe2gene.get(s), new ArrayList<String>());
			}
				gene2probe.get(probe2gene.get(s)).add(s);
			}
		Set <String> tmp3 = new HashSet<String>(gene2probe.keySet());
		List <String> tmp4 = new ArrayList<String>();
		for(String s: tmp3){
			if(gene2probe.get(s).size() <= 2){
				gene2probe.remove(s);	
			}else{
				tmp4.addAll(gene2probe.get(s));
			}	
		}
		R = R.getSubMatByRow(tmp4);
		rowSize = R.rowSize();
		nullDist = new HashMap<Integer, List<Double>>();
		setStatFunc2RankSum();	
		sampledIndex = new ArrayList<Integer>(); 
		for(int i = 1; i <= rowSize; i++){
			sampledIndex.add(i);
		}
	}
	
	private void getZscore(){
		zscore = true;
		nullDistSize = 1000;
	}
	
	static public MyMat rankMatrix(MyMat M, boolean ascending){
		MyMat R = new MyMat(M.getRowNames(), M.getColNames());
		for(String c: M.getColNames()){
			Map <String, Double> tmp = M.getColMap(c);
			List <String> sorted;
			if(ascending){
				sorted  = MyFunc.sortKeysByAscendingOrderOfValues(tmp);
			}else{
				sorted  = MyFunc.sortKeysByDescendingOrderOfValues(tmp);
			}
			for(int i=0;i<sorted.size();i++){
				R.set(sorted.get(i), c, i+1);
			}
		}
		return R;
	}
	
	private void rankMatrix(){
		S = rankMatrix(R, ascending);
	}
	
	// the smaller, the more significant 
	static public double rankProduct(List <Integer> L){
		double d = 1;	
		for(Integer s: L){
			d *= s;
		}
		return Math.pow(d, -L.size());
	}
	
	// the smaller, the more significant
	static public double rankProductForTopK(List <Integer> L, int k){
		double d = 1;	
		Collections.sort(L);
		for(int i = 0; i < k; i++){
			d *= L.get(i);
		}
		return Math.pow(d, -k);
	}
	
	// the smaller, the more significant 
	static public double rankSum(List <Integer> L){
		double d = 0;	
		for(Integer s: L){
			d += s;
		}
		return d/L.size();
	}
	
	// the smaller, the more significant
	static public double rankSumForTopK(List <Integer> L, int k){
		double d = 0;	
		Collections.sort(L);
		for(int i = 0; i < k; i++){
			d += L.get(i);
		}
		return d/k;
	}
	
	//the greater, the more significant
	static public double ksScore(List <Integer> L, int n){
		Set <Integer> S = new HashSet<Integer>(L);
		double d = 0;
		int k1 = L.size();
		int k2 = n-L.size();
		double  max = 0;	
		for(int i = 1; i <= n;i++){
			if(S.contains(i)){
				d +=  k2;
			}else{
				d -= k1;
			}
			if(d > max){
				max = d;
			}
		}
		return  max;
	}
	
	interface StatFunc {
		double get(List <Integer> L);
	}		
	
	StatFunc statFunc;
	StatFuncType statFuncType;
	public static class StatFuncType  implements Serializable {
		private static final long serialVersionUID = -1140169245215004688L;
		private final String name;
		private StatFuncType(String name){this.name = name;}
		public StatFuncType(StatFuncType b){this.name = b.name;}
		@Override
		public String toString(){return name;}
		public boolean equals(StatFuncType b){return name.equals(b.toString());};
		public static final StatFuncType RANKSUM = new StatFuncType("ranksum");
		public static final StatFuncType RANKSUM2 = new StatFuncType("ranksum2");
		public static final StatFuncType RANKPRODUCT = new StatFuncType("rankproduct");
		public static final StatFuncType RANKPRODUCT2= new StatFuncType("rankproduct2");
		public static final StatFuncType KSSCORE = new StatFuncType("ksscore");
	}
			
	private void setStatFunc2RankSum(){
		statFuncType = StatFuncType.RANKSUM;
		statFunc =new StatFunc(){
			public double get(List<Integer> L){
				return  -rankSum(L);
			}
		};
	}
	
	private void setStatFunc2RankSum2(){
		statFuncType = StatFuncType.RANKSUM2;
		statFunc =new StatFunc(){
			public double get(List<Integer> L){
				return  -rankSumForTopK(L,2);
			}
		};
	}
	
	private void setStatFunc2RankProduct(){
		statFuncType = StatFuncType.RANKPRODUCT;
		statFunc =new StatFunc(){
			public double get(List<Integer> L){
				return  -rankProduct(L);
			}
		};
	}
	
	private void setStatFunc2RankProduct2(){
		statFuncType = StatFuncType.RANKPRODUCT2;
		statFunc =new StatFunc(){
			public double get(List<Integer> L){
				return  -rankProductForTopK(L,2);
			}
		};
	}
	
	private void setStatFunc2ksScore(){
		statFuncType = StatFuncType.KSSCORE;
		statFunc =new StatFunc(){
			public double get(List<Integer> L){
				return  ksScore(L,rowSize);
			}
		};
	}
	
	private void calculateStatistics(){
		T = new MyMat(new ArrayList<String>(gene2probe.keySet()),R.getColNames());
		for(String s: T.getColNames()){
			for(String g: T.getRowNames()){
				List <Integer> L = new ArrayList<Integer>();
				for(String p: gene2probe.get(g)){
					L.add((int)S.get(p, s));
				}
				T.set(g, s, statFunc.get(L));
			}
		}
	}
	
	private void calculatePvalues(){
		P = new MyMat(new ArrayList<String>(gene2probe.keySet()),R.getColNames());	
		for(String s: T.getColNames()){
			for(String g: T.getRowNames()){
				double stat = T.get(g, s);
				int k = gene2probe.get(g).size();
				List <Double> nullDist = getNullDist(k);
				double p = 0;
				for(Double d: nullDist){
					if(d > stat){
						p++;
					}
				}
				if(p ==0){
					p = 1;
				}
				p /= nullDist.size();
				p = -Math.log10(p);
				P.set(g,s,p);
			}
		}
	}
	
	private void calculateZscore(){
		P = new MyMat(new ArrayList<String>(gene2probe.keySet()),R.getColNames());	
		for(String s: T.getColNames()){
			for(String g: T.getRowNames()){
				double stat = T.get(g, s);
				int k = gene2probe.get(g).size();
				List <Double> nullDist = getNullDist(k);
				
				double m = MyFunc.mean(nullDist);
				double sd = MyFunc.sd(nullDist);
				double z = (stat-m)/sd;
				P.set(g,s,z);
			}
		}
	}
	
	
	private List<Double>getNullDist(int k){
		if(nullDist.containsKey(k)){
			return nullDist.get(k);
		}else{ 
			List <Double> tmp = new ArrayList<Double>();
			while(tmp.size()< nullDistSize){
				tmp.add(statFunc.get(MyFunc.sampleForSmallN(sampledIndex,k)));
			}
			nullDist.put(k, tmp);
			return tmp;
		}
	}
	
	
	public void perform(){
		System.err.println("rank matrix....");
		rankMatrix();
		System.err.println("calculate " + statFuncType + " statistics....");
		calculateStatistics();
		if(!zscore){
			System.err.println("calculate pvalues....");
			calculatePvalues();
		}else{
			System.err.println("calculate zscore....");
			calculateZscore();
		}
	}
	
	                   
	public MyMat getPvalueMatrix(){
		return P;
	}
	private void print(){
		System.out.print(P);
	}
	
	private void print(String outfile) throws IOException{
		P.print(outfile);
	}
	
	
	private void setStatFuncType(Integer I){
		switch(I){
			case 1:
				setStatFunc2RankProduct();
			case 2:
				setStatFunc2RankProduct2();
			case 3:
				setStatFunc2RankSum();
			case 4:	
				setStatFunc2RankSum2();
			case 5:
				setStatFunc2ksScore();
			default:				
		}
	}
	
	
	
	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("s", "stat", true, "stat func type (1:rankproduct, 2:rankproduct, 3:ranksum, 4:ranksum, 5:ksscore)");
		options.addOption("n", "ndsize", true, "null dist size");
		options.addOption("d", "dec", false, "rank for the descending order");
		options.addOption("o", "outfile", true, "outfile");
		options.addOption("z", "zscore", false, "get z score");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabFile chipFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] tabFile chipFile", options);
			return;
		}
	
		CollapseMatrixRow CMR = new CollapseMatrixRow(new MyMat(argList.get(0)), MyFunc.readStringStringMap(argList.get(1)));
		if(commandLine.hasOption("s")){
			CMR.setStatFuncType(Integer.valueOf(commandLine.getOptionValue("s")));
		}
		if(commandLine.hasOption("n")){
			CMR.nullDistSize = Integer.valueOf(commandLine.getOptionValue("n"));
		}
		if(commandLine.hasOption("d")){
			CMR.ascending = false;
		}
		if(commandLine.hasOption("z")){
			CMR.getZscore();
		}
		CMR.perform();
		if(commandLine.hasOption("o")){
			CMR.print(commandLine.getOptionValue("o"));
		}else{
			CMR.print();
		}
	}
	
	
}
