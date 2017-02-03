package utility;
import java.util.*;
import java.io.*;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import sun.reflect.Reflection;
import eem.ExpressionModuleSet;


public class GeneSetOverlapTest{
	private Map <String, List <String>> geneSetMap1;
	private Map <String, List <String>> geneSetMap2;
	
	private List <String> bgGenes;
	
	private int bgGeneNumber = 10000;
	
	private MyMat Pvalue;
	private MyMat Qvalue;
	
	private Double PvalueCutoff = null;
	private Double QvalueCutoff = null;	
	private boolean minusLogScale = false;	
	
	
	public GeneSetOverlapTest(Map <String, List <String> > geneSetMap1, Map <String, List <String> > geneSetMap2) {
		bgGenes = new ArrayList<String>();
		setGeneSetMap1(geneSetMap1);		setGeneSetMap2(geneSetMap2);
		
	}	
	

	public GeneSetOverlapTest(Map <String, List <String> > geneSetMap1) {
		bgGenes = new ArrayList<String>();
		setGeneSetMap1(geneSetMap1);	
	}	
	
	public void setBgGenes(List <String> Gene){
		bgGenes = new ArrayList<String>(Gene);
		bgGeneNumber = bgGenes.size();	
		if(geneSetMap1 != null && !geneSetMap1.isEmpty()){
			Map <String, List<String>> tmp = new HashMap<String, List<String>>();	
			for(Map.Entry<String, List<String>> e: geneSetMap1.entrySet()){
				List <String> tmp2 = MyFunc.isect(bgGenes, e.getValue());
				if(!tmp2.isEmpty()){
					tmp.put(e.getKey(),tmp2);
				}			
			}
			geneSetMap1 = tmp;
		}
		if(geneSetMap2 != null && !geneSetMap2.isEmpty()){
			Map <String, List<String>> tmp = new HashMap<String, List<String>>();	
			for(Map.Entry<String, List<String>> e: geneSetMap2.entrySet()){
				List <String> tmp2 = MyFunc.isect(bgGenes, e.getValue());
				if(!tmp2.isEmpty()){
					tmp.put(e.getKey(),tmp2);
				}			
			}
			geneSetMap2 = tmp;
		}
		
	}
	
	public void setBgGeneNumber(int n){
			bgGeneNumber = n;
	}
	
		
	public void setPvalueCutoff(double d){
		PvalueCutoff = d;	
	}
	public void setQvalueCutoff(double d){
		QvalueCutoff = d;	
	}
	public void outputInMinusLogScale(){
		minusLogScale = true;
	}
	
	
	public double getPvalue(String geneSetID1, String geneSetID2){
		return Pvalue.get(geneSetID1, geneSetID2);	
	}
	
	public double getQvalue(String geneSetID1, String geneSetID2){
		return Qvalue.get(geneSetID1, geneSetID2);	
	}
	
	public Map<String, Double> getPvalueMap(){
		return  asMap(Pvalue);
	}
	public Map<String, Double> getQvalueMap(){
		return asMap( Qvalue);
	}
	
	
	public void setGeneSetMap1(Map <String, List <String> > geneSetMap){
		if(bgGenes.isEmpty()){
			geneSetMap1 = geneSetMap;
		}else{
			Map <String, List<String>> tmp = new HashMap<String, List<String>>();	
			for(Map.Entry<String, List<String>> e: geneSetMap.entrySet()){
				List <String> tmp2 = MyFunc.isect(bgGenes, e.getValue());
				if(!tmp2.isEmpty()){
					tmp.put(e.getKey(),tmp2);
				}			
			}
			geneSetMap1 = tmp;
		}	
	}
	
	public void setGeneSetMap2(Map <String, List <String> > geneSetMap){
		if(bgGenes.isEmpty()){
			geneSetMap2 = geneSetMap;
		}else{
			Map <String, List<String>> tmp = new HashMap<String, List<String>>();	
			for(Map.Entry<String, List<String>> e: geneSetMap.entrySet()){
				List <String> tmp2 = MyFunc.isect(bgGenes, e.getValue());
				if(!tmp2.isEmpty()){
					tmp.put(e.getKey(),tmp2);
				}			
			}
			geneSetMap2 = tmp;
		}	
	}
	
	
	private  Map <String, Double> asMap(MyMat M){
		Map <String, Double> m = new HashMap<String, Double>();
		if(geneSetMap2 != null){
		for(String s: M.rowname){	
			for(String t: M.colname){
				m.put(s + "\t" + t,M.get(s,t));	
			}
		}
		}else{
			List <String> tmp = M.getRowNames();
			for(String s: tmp){
				for(String t: tmp){
					if(s.compareTo(t)> 0){
						m.put(s + "\t" + t,M.get(s,t));
					}		
				}
			}
		}
		return m;		
	}
	public void calculatePvalue(){
		if(geneSetMap2 != null ){
			Pvalue = new MyMat(new ArrayList<String>(geneSetMap1.keySet()), new ArrayList<String>(geneSetMap2.keySet()));
			for(Map.Entry<String, List<	String>> e1: geneSetMap1.entrySet()){
				for(Map.Entry<String, List<String>> e2: geneSetMap2.entrySet()){
					int isect = (MyFunc.isect(e1.getValue(), e2.getValue())).size();
					double P = (isect == 0) ? 1 : MyFunc.calculatePvalueForSetOverlap(bgGeneNumber, e1.getValue().size(), e2.getValue().size(),isect );
					System.err.println( e1.getKey() + "(" +  e1.getValue().size() + ")" + " vs " +  e2.getKey() + "(" +  e2.getValue().size() + ")" +  ":\t" + "isect=" + isect + "\t" + "p=" + P );
					Pvalue.set(e1.getKey(), e2.getKey(), P);
				}
			}
		}else{
			Pvalue = new MyMat(new ArrayList<String>(geneSetMap1.keySet()), new ArrayList<String>(geneSetMap1.keySet()));
			List <String> tmp = new ArrayList<String>(geneSetMap1.keySet());
			int isect;
			double P;
			for(int i = 0, n = tmp.size(); i <n ;i++){
				Pvalue.set(tmp.get(i), tmp.get(i), 1);		
				for(int j = 0; j <= i ;j++){
					isect = (MyFunc.isect(geneSetMap1.get(tmp.get(i)), geneSetMap1.get(tmp.get(j))).size());
					P = (isect == 0) ? 1 : MyFunc.calculatePvalueForSetOverlap(bgGeneNumber, geneSetMap1.get(tmp.get(i)).size(), geneSetMap1.get(tmp.get(j)).size(),isect );
					System.err.println( tmp.get(i) + "(" +  geneSetMap1.get(tmp.get(i)).size() + ")" + " vs " +  tmp.get(j) + "(" +  geneSetMap1.get(tmp.get(j)).size() + ")" +  ":\t" + "isect=" + isect + "\t" + "p=" + P );
					Pvalue.set(tmp.get(i), tmp.get(j), P);
					Pvalue.set(tmp.get(j), tmp.get(i), P);
				}
			}
		}

		
	}
	
	
	
	public void calculateQvalue(){
		if(geneSetMap2 != null ){
			Qvalue = new MyMat(new ArrayList<String>(geneSetMap1.keySet()), new ArrayList<String>(geneSetMap2.keySet()));
		}else{
			Qvalue = new MyMat(new ArrayList<String>(geneSetMap1.keySet()), new ArrayList<String>(geneSetMap1.keySet()));
		}
		Map <String, Double> PvalueMap = Pvalue.asMap();
		Map <String, Double> QvalueMap = MyFunc.calculateStoreyQvalue(PvalueMap);
		for(Map.Entry<String, Double> e: QvalueMap.entrySet()){
			List <String> tmp = MyFunc.split("\t", e.getKey());
			Qvalue.set(tmp.get(0), tmp.get(1), e.getValue());
			if(geneSetMap2 == null ){
				Qvalue.set(tmp.get(1), tmp.get(0), e.getValue());
			}	
		}
		if(geneSetMap2 == null ){
			List <String> tmp = new ArrayList<String>(geneSetMap1.keySet());
			for(int i = 0, n = tmp.size(); i <n ;i++){
				Qvalue.set(tmp.get(i), tmp.get(i), 0);	
			}				
		}
	}
	
	
	public String toString(){
		StringBuffer S = new StringBuffer("\t\tP value\tQ value\n");
		Map <String, Double> PvalueMap = getPvalueMap();
		Map <String, Double> QvalueMap = getQvalueMap();
		List<String> keys =  MyFunc.sortKeysByAscendingOrderOfValues(PvalueMap);
		for(String s: keys){
			List <String>  tmp = new ArrayList<String>();
			tmp.add(s);
			double p; 
			if(minusLogScale){
				p = (PvalueMap.get(s)==0)?Double.MAX_VALUE: -Math.log10(PvalueMap.get(s));
			}else{
				p = PvalueMap.get(s);
			}
			if(PvalueCutoff != null && (minusLogScale?(p < PvalueCutoff):(p > PvalueCutoff))){
				continue;
			}
			double q; 
			if(minusLogScale){
				q = (QvalueMap.get(s)==0)?Double.MAX_VALUE: -Math.log10(QvalueMap.get(s));
			}else{
				q = QvalueMap.get(s);
			}
			if(QvalueCutoff != null &&  (minusLogScale?(q < QvalueCutoff):(q > QvalueCutoff))){
				continue;
			}
			tmp.add(Double.toString(p));
			tmp.add(Double.toString(q));
			S.append(MyFunc.join("\t", tmp) + "\n");
		}
		return S.toString();
	}
	
	public MyMat getMinusLogPvalueMatrix(){	
		double max = 20;	
			MyMat tmp = new MyMat(Pvalue.getRowNames(), Pvalue.getColNames());
			for(String s: Pvalue.getRowNames()){
				for(String t: Pvalue.getColNames()){
					double p = (Pvalue.get(s, t)==0)?Double.MAX_VALUE:-Math.log10(Pvalue.get(s, t));
					if(p>max){
						p=max;
					}
					tmp.set(s, t, p);
				}
			}
			return tmp;
	}
	
	
	public static Map<String, List<String>> setGeneSetFromList(String infile) throws IOException{
		Map<String, List<String>> geneSet = new HashMap<String, List<String>>();
		geneSet.put(infile, MyFunc.readStringList2(infile));
		return geneSet;
	}
	
	public static void main(String[] args) throws Exception {
		Options options = new Options();
		options.addOption("L", "list", true, " first input from list file");
		options.addOption("w", "wlist", false, "both from list file");
		options.addOption("b", "bg", true, "# of background genes");
		options.addOption("p", "pcutoff", true, "overlap p-value cutoff");
		options.addOption("q", "qcutoff", true, "overlap q-value cutoff");
		options.addOption("P", "pmat", false, "get minus log p-value matrix");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		
		
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] gmt_file gmt_file", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2 && argList.size() != 1){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] gmt_file gmt_file", options);
			return;
		}
		
		GeneSetOverlapTest GSO = null;
		if(!(commandLine.hasOption("L") || commandLine.hasOption("w"))){
			if(argList.size() != 1){
				GSO = new GeneSetOverlapTest(MyFunc.readGeneSetFromGmtFile(argList.get(0)), MyFunc.readGeneSetFromGmtFile(argList.get(1)));
			}else{
				GSO = new GeneSetOverlapTest(MyFunc.readGeneSetFromGmtFile(argList.get(0)));
			}
		}else{
			if(commandLine.hasOption("L") & argList.size() == 1){
				GSO = new GeneSetOverlapTest(setGeneSetFromList(commandLine.getOptionValue("L")), MyFunc.readGeneSetFromGmtFile(argList.get(0)));
			}else if (commandLine.hasOption("w") & argList.size() == 2){
				GSO = new GeneSetOverlapTest(setGeneSetFromList(argList.get(0)), setGeneSetFromList(argList.get(1)));
			}else{
				formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] gmt_file gmt_file", options);
				return;
			}	
		}
		if(commandLine.hasOption("b")){
			GSO.setBgGeneNumber(Integer.valueOf(commandLine.getOptionValue("b")));
		}
		if(commandLine.hasOption("p")){
			GSO.setPvalueCutoff(Double.valueOf(commandLine.getOptionValue("p")));
		}
		if(commandLine.hasOption("q")){
			GSO.setQvalueCutoff(Double.valueOf(commandLine.getOptionValue("q")));
		}
		
		GSO.calculatePvalue();
		if(commandLine.hasOption("P")){
			System.out.print(GSO.getMinusLogPvalueMatrix());
		}else{
			GSO.calculateQvalue();
			System.out.print(GSO);
		}
	}			
	
}
