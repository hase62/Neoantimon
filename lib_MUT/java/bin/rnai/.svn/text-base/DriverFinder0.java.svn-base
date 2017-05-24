package rnai;

import sun.reflect.Reflection;
import utility.*;

import java.io.IOException;
import java.util.*;
import org.apache.commons.cli.*;

public class DriverFinder0 {
	MyMat R;
	MyMat T;
	MyMat C;
	MyMat S;
	Map <String, List<String>> gene2probe;	
	int minListSize = 6;

	double cutoff = 2; 
	
	
	public DriverFinder0(MyMat r, MyMat t, Map <String, String> probe2gene) throws Exception{
		
		if(MyFunc.isect(r.getColNames(), t.getColNames()).size() >= minListSize){
			R = r;
			T = t;
		}else if(MyFunc.isect(r.getRowNames(), t.getColNames()).size() >= minListSize){
			R = r;
			T = t;
			r.transpose();
		}else if(MyFunc.isect(r.getColNames(), t.getRowNames()).size() >= minListSize){
			R = r;
			T = t;
			t.transpose();
		}else if(MyFunc.isect(r.getRowNames(), t.getRowNames()).size() >= minListSize){
			R = r;
			T = t;
			r.transpose();
			t.transpose();
		}else{
			throw new Exception("ERR: no row and column ids");
		}	
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
	}
	
	
	private void getCorrelation(){
		MatrixCorrelation2 MC = null;
		try {
			MC = new MatrixCorrelation2(R,T);
		} catch (Exception e) {
			e.printStackTrace();
		}
		MC.calculatePvalue();
		C = MC.getMinusLogPMatrix();
	}
	
	private void convertProbe2gene(){
		S = new MyMat(new ArrayList<String>(gene2probe.keySet()), C.getColNames());
		for(String s: C.getColNames()){
			for(String g:gene2probe.keySet()){
				List <String> tmp2 = gene2probe.get(g);
				List <Double> tmp3 = new ArrayList<Double>();
				for(String p: tmp2){
					tmp3.add(C.get(p,s));	
				}
				
				Collections.sort(tmp3);
				double tmp4 = tmp3.get(1);
				Collections.reverse(tmp3);
				if(Math.abs(tmp4) < Math.abs(tmp3.get(1))){
					tmp4 = tmp3.get(1);
				}
				S.set(g, s,tmp4);
			}	
		}
	}
	
	
	private void filterOutputRow(double cutoff){
		S = filterMatrixRow(S,cutoff);
	}
	
	private void filterOutputColumn(double cutoff){
		S = filterMatrixColumn(S,cutoff);
	}
	
	
	static private MyMat filterMatrixRow(MyMat S, double cutoff){
		List <String> filtered = new ArrayList<String>(); 
		L:for(String s:S.getRowNames()){
			for(String t: S.getColNames()){
				if(Math.abs(S.get(s, t)) > cutoff){
					filtered.add(s);
					continue L;
				}
			}
		}
		if(filtered.size()>0){
			S = S.getSubMatByRow(filtered);
		}	
		return S;
	}
	
	static private MyMat filterMatrixColumn(MyMat S, double cutoff){
		List <String> filtered = new ArrayList<String>(); 
		L:for(String s:S.getColNames()){
			for(String t: S.getRowNames()){
				if(Math.abs(S.get(t, s)) > cutoff){
					filtered.add(s);
					continue L;
				}
			}
		}
		if(filtered.size()>0){
			S = S.getSubMatByRow(filtered);
		}
		return S;
	}
	
	private void print(){
		System.out.print(S);
	}
	
	private void print(String outfile) throws IOException{
		S.print(outfile);
	}
	
	
	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("r", "row", true, "filter row");
		options.addOption("c", "column", true, "filter column");
		options.addOption("o", "outfile", true, "outfile");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] rnaiTabFle targetTabFile chipFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 3 && argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] rnaiTabFle targetTabFile chipFile", options);
			return;
		}
		
		DriverFinder0 DF = new DriverFinder0(new MyMat(argList.get(0)),	new MyMat(argList.get(1)), MyFunc.readStringStringMap(argList.get(2)));
		DF.getCorrelation();
		DF.convertProbe2gene();
		if(commandLine.hasOption("c")){
			DF.filterOutputColumn( Double.valueOf(commandLine.getOptionValue("c")));
		}
		if(commandLine.hasOption("r")){
			DF.filterOutputRow(Double.valueOf(commandLine.getOptionValue("r")));
		}
		if(commandLine.hasOption("o")){
			DF.print(commandLine.getOptionValue("o"));
		}else{
			DF.print();
		}
	}
}
