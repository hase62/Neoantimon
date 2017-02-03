package snp;


import java.util.*;
import java.io.*;
import org.apache.commons.cli.*;

import utility.PoibinRNA;

import sun.reflect.Reflection;
import utility.*;


public class PACT {
	MyMat A;
	MyMat B;
	MyMat Pr;
	int m;
	int n;
	Dist C;
	Dist J; // Jaccard Index
	
	double jcutoff = 0.1;
	
	double posWeight = 1;
	double sampWeight = 1;
	
	boolean  exclusive = false;
	
	
	
	
	Map<String, Double> P;
	Map<String, Double> Q;
	List <Double> pp;
	double cutoff = 0.5;
	boolean countLess = false ;
	
	public PACT (MyMat M){
		A = M;
		n = A.rowSize();
		m = A.colSize();
	}
	
	
	
	
	public  void setCutoff(double d){
		cutoff = d;
	}
	public  void setCutoffByPercentile(double d){
		if(d > 1){
			d /= 100;
		}
		cutoff = MyFunc.percentile(A.asList(), d);
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
	
	private void binarizeMatrix(){
		B = new MyMat(A.getRowNames(), A.getColNames());
		if(!countLess){
			for(int i = 0; i < n; i++){
				for(int j = 0; j< m; j++){
					if(A.get(i, j) >= cutoff){
						B.set(i, j, 1.0);
					}
				}
			}
		}else{
			for(int i = 0; i < n; i++){
				for(int j = 0; j< m; j++){
					if(A.get(i, j) <= cutoff){
						B.set(i, j, 1.0);
					}
				}
			}	
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
	
	
	private void calculateParameters(){
		int a;
		List <Double> q;
		List <Double> r;
		a = (int)MyFunc.sum(B.asList());
		q = new ArrayList<Double>();
		r = new ArrayList<Double>();
		
		for(int i = 0; i < n; i++){
			q.add(MyFunc.sum(B.getRow(i))/a);
		}
		
		if(posWeight  != 1){
			q = reWeight(q,posWeight);
		}	
		
		
		for(int k = 0; k < m; k++){
			r.add(MyFunc.sum(B.getCol(k))/a);
		}
		if(sampWeight != 1){
			r = reWeight(r,sampWeight);
		}
		Pr = new MyMat(A.getRowNames(), A.getColNames());
		for(int i = 0; i < n; i++){
			for(int k = 0;k < m; k++){
				double p = a*q.get(i)*r.get(k);
				 if(p > 1){
					System.err.println(p); 
					//p = 0.99999;
				 }
				Pr.set(i, k, p);
			}
		}
		
	}
	
	private List <Double> getPoiBinParameters(int i, int j){
		List <Double> tmp = new ArrayList <Double>();
		if(!exclusive){
			for(int k = 0; k < m; k++){
				double p  = Pr.get(i,k)*Pr.get(j,k);
				tmp.add(p);
			}
		}else{
			for(int k = 0; k < m; k++){
				double p  = Pr.get(i,k)*(1-Pr.get(j,k)) + (1-Pr.get(i,k))*Pr.get(j,k);;
				tmp.add(p);
			}
		}
		return tmp;
	}
	
	
	
	private void getCoaberrationCount(){
		C = new Dist(A.getRowNames());
		J = new Dist(A.getRowNames());
		if(!exclusive){
			for(int i=0; i<n;i++){
				for(int j=0; j<i;j++){
					double and = 0;
					double or = 0;
					for(int k=0; k<m;k++){
						if(B.get(i,k) == 1 & B.get(j,k) == 1){
							and++;
						}
						if(B.get(i,k) == 1 | B.get(j,k) == 1){
							or++;
						}
					}
					C.set(i, j, and);
					J.set(i, j, (or==1)?-1:and/or);
				
				}
			}
		}else{
			for(int i=0; i<n;i++){
				for(int j=0; j<i;j++){
					double or = 0;
					double exclusiveOr  = 0;
					int i2 = 0;
					int j2 = 0;
					for(int k=0; k<m;k++){
						if(B.get(i,k) == 1 | B.get(j,k) == 1){
							or++;
						}
						if(Math.abs(B.get(i,k)-B.get(j,k))==1){
							exclusiveOr++;
						}
						if(B.get(i,k) == 1){
							i2++;
						}
						if(B.get(j,k) == 1){
							j2++;
						}
					}
					C.set(i, j, exclusiveOr);
					J.set(i, j, (i2==0 | j2==0)?-1:exclusiveOr/or);
				}
			}
			
		}
	}
	
	
	private  void calculatePvalues(){
		P =  new  HashMap<String, Double>();
		PoibinRNA PB = new PoibinRNA();
		for(int i=0; i< C.size();i++){
			for(int j=0; j<i;j++){
				double pvalue;
				if(J.get(i, j) > jcutoff ){
					int k = (int)C.get(i,j);
					List <Double> p = getPoiBinParameters(i, j);
					pvalue = 1- PB.getCdf(k,p);
				}else{
					pvalue=1;
				}
				
				P.put(C.getNames().get(i) + "\t"+  C.getNames().get(j), pvalue);
			}
		}
	}
		
	private  void calculateQvalues(){
		Q = MyFunc.calculateQvalue(P);
		for(String p: Q.keySet()){ 
			Q.put(p, Q.get(p)>1?1:Q.get(p));
		}
	}
	
	private Dist getMinusLogPmatrix(){
		Map<String, Double> p = getMinusLogPvalues();
		Dist P = new Dist (C.getNames());
		for(int i = 0; i < C.size(); i++){
			for(int j = 0; j < i; j++){
				P.set(i, j, p.get(C.getNames().get(i) + "\t" + C.getNames().get(j)));
			}
		}
		return P;
	}
	
	public void perform(){
		System.err.println("binarize  matrix....");
		binarizeMatrix();
		
		System.err.println("calculate parameters....");
		calculateParameters();
		
		System.err.println("count coaberrations....");
		getCoaberrationCount();
		
		System.err.println("calculate pvalues....");
		calculatePvalues();
		
		System.err.println("calculate qvalues....");
		calculateQvalues();
	}
	
	public Map <String, Double> getMinusLogQvalues(){
		Map <String, Double> minusLogQvalues = new LinkedHashMap<String, Double>();
		double  tmp = 1;
		for(String p: Q.keySet()){
			if(Q.get(p) != 0 &  Q.get(p) < tmp){
				tmp = P.get(p);
			}
		}
		double max = - Math.log10(tmp);
		
		for(String p: Q.keySet()){
			tmp = Q.get(p);
			if(tmp==1){
				tmp=0;
			}else if(tmp==0){
				tmp = max;
			}else{
				tmp = - Math.log10(tmp);
			}
			minusLogQvalues.put(p, tmp);
		}
		
		return minusLogQvalues;
	}
	
	public Map <String, Double> getMinusLogPvalues(){
		Map <String, Double> minusLogPvalues = new LinkedHashMap<String, Double>();
		
		double  tmp = 1;
		for(String p: P.keySet()){
			if(P.get(p) != 0 &  P.get(p) < tmp){
				tmp = P.get(p);
			}
		}
		double max = - Math.log10(tmp);
		
		for(String p: P.keySet()){
			tmp = P.get(p);
			if(tmp==1){
				tmp=0;
			}else if(tmp==0){
				tmp = max;
			}else{
				tmp = - Math.log10(tmp);
			}
			minusLogPvalues.put(p, tmp);
		}
		return minusLogPvalues;
	}
	
	public static void main(String [] args) throws Exception{
		Options options = new Options();
		options.addOption("e", "exc", false, "test exclusive or");
		options.addOption("m", "tab", false, "read tab file");
		options.addOption("r", "r", true, "print parameter matrix");
		options.addOption("c", "cutoff", true, "cutoff for aberration");
		options.addOption("o", "outfile", true, "output file name");
		options.addOption("t", "thinpi", true, "thin down probe info");
		options.addOption("n", "nump", true, "# of psuedo probes");
		options.addOption("l", "less", false, "count less than cutoff");
		options.addOption("p", "pcutoff", true, "pvalue cutoff");
		options.addOption("q", "qcutoff", true, "qvalue cutoff");
		options.addOption("d", "diffchrom", false, "different chromosomes");
		options.addOption("C", "cutbyper", true, "cutoff by percentile");
		options.addOption("b", "bed", false, "gene coordinate file from bed file");
		options.addOption("T", "tab", false, "print tab");
		options.addOption("P", "suppos", true, "suppress position-wise weight");
		options.addOption("S", "samppos", true, "suppress sample-wise weight");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine = null;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile  probeTsvFile", options);
			System.exit(1);
		}
		List <String> argList = commandLine.getArgList();
		if(!(argList.size() == 1 | argList.size() == 2)){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] SegFile probeTsvFile", options);
			System.exit(1);
		}
		
		ProbeInfo PI = null;
		PACT  PACT;
		if(commandLine.hasOption("m")){
			PACT = new PACT(new MyMat(argList.get(0)));
		}else{
			SegmentContainerMap SCM = new SegmentContainerMap(argList.get(0));
			List<String> sample = new ArrayList<String>(SCM.keySet());
			if(argList.size() == 2 ){
				if(!commandLine.hasOption("b")){
					PI = 	ProbeInfo.getProbeInfoFromTsvFile(argList.get(1));
				}else{
					PI = GeneInfo.getGeneInfoFromBedFile(argList.get(1)).toProbeInfo();
				}
				if(commandLine.hasOption("t")){
					PI.thinDown(Integer.valueOf(commandLine.getOptionValue("t")));
				}
				PI.filter(SCM.get(sample.get(0)));
			
		}else{
			int n = 100;
			if(commandLine.hasOption("n")){
				n = Integer.valueOf(commandLine.getOptionValue("n"));
			}	
			PI = SCM.generatePsuedoProbeInfo(n);
		}
		PACT = new PACT(SCM.toMyMat(PI));
		}
		if(commandLine.hasOption("S")){
			PACT.setSampWeight(Double.valueOf(commandLine.getOptionValue("S")));
		}
		if(commandLine.hasOption("P")){
			PACT.setPosWeight(Double.valueOf(commandLine.getOptionValue("P")));
		}
		if(commandLine.hasOption("l")){
			PACT.countLess = true;
		}
		
		if(commandLine.hasOption("c")){
			PACT.setCutoff(Double.valueOf(commandLine.getOptionValue("c")));
		}
		if(commandLine.hasOption("l")){
			PACT.countLess = true;
		}
		if(commandLine.hasOption("C")){
			PACT.setCutoffByPercentile(Double.valueOf(commandLine.getOptionValue("C")));
		}
		if(commandLine.hasOption("e")){
			PACT.exclusive = true;
		}
			
		PACT.perform();		

		if(commandLine.hasOption("r")){
			PACT.Pr.print(commandLine.getOptionValue("r"));
		}
		
		Writer os;
		if(commandLine.hasOption("o")){
			 os = new BufferedWriter(new FileWriter(commandLine.getOptionValue("o")));	 
		}else{
			os = new PrintWriter(System.out);
		}
		
		
		if(commandLine.hasOption("T")){
			os.write(PACT.getMinusLogPmatrix().toStringInTabFormat());
			os.flush();
			System.exit(0);
		}
		
		double pcutoff = -1;
		double qcutoff = -1;
		boolean diffchrom = false;
		
		
		if(commandLine.hasOption("p")){
			pcutoff = Double.valueOf(commandLine.getOptionValue("p"));
		}
		if(commandLine.hasOption("q")){
			qcutoff = Double.valueOf(commandLine.getOptionValue("q"));
		}
		
		if(commandLine.hasOption("d")){
			if(PI != null){
				diffchrom = true;
			}
		}
		
		os.write("Probe1" + "\t" +  "Probe2"  + "\t" + "JaccardIndex" + "\t" + "pvalue" + "\t" + "qvalue" + "\n");
		Map<String, Double> p = PACT.getMinusLogPvalues();
		Map<String, Double> q = PACT.getMinusLogQvalues();
			for(String s: MyFunc.sortKeysByDescendingOrderOfValues(p)){
				if(p.get(s) > pcutoff & q.get(s) > qcutoff){
					List<String> tmp = MyFunc.split("\t", s);
					if(diffchrom ){
						if(PI.chr(tmp.get(0)) ==  PI.chr(tmp.get(1))){
							continue;
						}
					}
					os.write(s  +  "\t" +  PACT.J.get(tmp.get(0), tmp.get(1))  + "\t"+ p.get(s) + "\t" + q.get(s) + "\t" + "\n");
				}
			}
			os.flush();
	}

	
	
}
