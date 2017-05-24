package mutation;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.DataFormatException;

import mutation.HotNeighborTest.geneComp;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.math3.distribution.BinomialDistribution;

import sun.reflect.Reflection;
import utility.*;

public class HotGenesetTest{
	private Map <String,Integer> mutationCount;
	private Map <String,Integer> mutationSumInGeneset;
	private Map <String,Integer> genesetSize;
	private Map <String, List<String>> geneset;
	private Map <String, List<String>> hotGenesetMember;
	private Map <String,Double> pvalue;
	private int totalMutationCount;
	private int geneUniverseSize = 20000;
	private double pcutoff = 3;
	private double pmax = 20;
	private Map <String,Double> qvalue;
	
	public HotGenesetTest(){}
	
	public void setPcutoff(double d){
		pcutoff = d;
	}
	public void applyCeiling2MutationCount(int i){
		for(String s: mutationCount.keySet()){
			if(mutationCount.get(s) > i){
				mutationCount.put(s,i);
			}
		}
		setTotalMutationCount();
	}
	
	public void cutSmallMutationCount(int i){
		Set <String> tmp = new HashSet<String>(mutationCount.keySet());
		for(String s: tmp){
			if(mutationCount.get(s) < i){
				mutationCount.remove(s);
			}
		}
		setTotalMutationCount();
	}
	
	private void setTotalMutationCount(){
		totalMutationCount = 0;
		for(String k: mutationCount.keySet()){
			totalMutationCount += Integer.valueOf(mutationCount.get(k));
		}
	}
	
	public void setMutationCount(String infile) throws IOException, DataFormatException{
		Map<String, String> tmp = MyFunc.readStringStringMap(infile);
		mutationCount = new HashMap<String, Integer>();
		for(String k: tmp.keySet()){
			mutationCount.put(k, Integer.valueOf(tmp.get(k)));
		}
		setTotalMutationCount();
	}
	
	public void setMutationCount(Map<String, Integer> count){
		mutationCount = new HashMap<String, Integer>(count);	
		setTotalMutationCount();
	}
	
	public void setGeneSet(String infile) throws IOException, DataFormatException{
		geneset =  MyFunc.readGeneSetFromGmtFile(infile);
	}
		
	public void setGeneSet(Map <String, List<String>> geneset){
		this.geneset =  new HashMap <String, List<String>>();
	}
	
	private void getMutationSumInGeneset(){
		mutationSumInGeneset = new  HashMap <String,Integer>();
		genesetSize = new  HashMap <String,Integer>();
		pvalue =  new HashMap <String,Double>();
		hotGenesetMember = new HashMap <String, List<String>>();
		for(String gs: geneset.keySet()){
			List <String> member = geneset.get(gs);
			genesetSize.put(gs, member.size());
			for(String g: member){
				if(mutationCount.containsKey(g)){
					if(!mutationSumInGeneset.containsKey(gs)){
						mutationSumInGeneset.put(gs, mutationCount.get(g)); 
						List <String> tmp = new ArrayList<String>();
						tmp.add(g);
						hotGenesetMember.put(gs, tmp);
					}else{
						mutationSumInGeneset.put(gs, mutationSumInGeneset.get(gs) + mutationCount.get(g));
						hotGenesetMember.get(gs).add(g);
					}
				}
			}
		}
	}
	
	private void calculatePvalues(){
		pvalue = new HashMap<String, Double>();
		for(String gs: geneset.keySet()){
			pvalue.put(gs, calculatePvalue(gs));
		}
	}
	
	private  double calculatePvalue(String gs){
		int n =  totalMutationCount;
		if(!mutationSumInGeneset.containsKey(gs)){
			return 1;
		}
		int k = mutationSumInGeneset.get(gs);
		double pvalue = 1;
		if(n>0 & k > 0){
			double p = (double)(genesetSize.get(gs))/geneUniverseSize;
			BinomialDistribution BN = new BinomialDistribution(n, p); 
			try {
				pvalue = 1- BN.cumulativeProbability(k); 
			}catch (Exception e){
				e.printStackTrace(); 
			}
		}
		pvalue = -Math.log10(pvalue);
		if(pvalue > pmax){
			pvalue = pmax;
		}
		return pvalue;
	}
	
	public void printResult(){
		System.out.println("gene\tpvalue\tqvalue\tmutationSumInGeneset\tgenesetSize\tmutatedMember");
		for(String gs: MyFunc.sortKeysByDescendingOrderOfValues(pvalue)){
			if(pvalue.get(gs) >= pcutoff){
				System.out.print(gs + "\t"  + pvalue.get(gs) + "\t"+ qvalue.get(gs) + "\t"+  mutationSumInGeneset.get(gs) + "\t" + genesetSize.get(gs)  + "\t");
				List <String> tmp = hotGenesetMember.get(gs);
				java.util.Collections.sort(tmp, new geneComp());
				List <String> tmp2 = new ArrayList<String>(); 
				for(String hn: tmp){
					tmp2.add(hn + "[" +  mutationCount.get(hn) + "]");
				}
				System.out.println(tmp2.isEmpty()?"":MyFunc.join(" ", tmp2));
			}
		}
	}
	
	public void printResultSimple(){
		for(String gene: MyFunc.sortKeysByDescendingOrderOfValues(pvalue)){
			if(pvalue.get(gene) >= pcutoff){
				System.out.println(gene  + "\t"  + pvalue.get(gene) + "\t"+qvalue.get(gene)); 
			}
		}
	}
	
	public class geneComp implements Comparator<String> {
        public int compare(String g1, String g2) {
        	if(mutationCount.get(g1) < mutationCount.get(g2)){
        		return 1;
        	}else if (mutationCount.get(g1) > mutationCount.get(g2)){
        		return -1;
        	}else{
				return 0;
        	}
        }
	}
	public void perform(){
		getMutationSumInGeneset();
		calculatePvalues();
		calculateQvalues();
	}
	
	private void calculateQvalues(){
		List <String> keys = MyFunc.sortKeysByDescendingOrderOfValues(pvalue);
		int n = pvalue.size();
		int i ;
		qvalue = new HashMap<String, Double>();
		for(i = 0; i < keys.size(); i++){
			double q = (n * Math.pow(10, -pvalue.get(keys.get(i)))) / (i+1) ; 
			if(q>1){ 
				q=1;
			}
			q = -Math.log10(q);
			if(q > pmax){
				q = pmax;
			}
			qvalue.put(keys.get(i), q);
	   	}	
	}

	public static void main(String [] args) throws Exception{	
		Options options = new Options();
		options.addOption("p", "pcut", true,  "pvalue cutoff");
		options.addOption("s", "simple", false,  "print simple result");
		options.addOption("c", "ceiling", true,  "mutation count ceiling");
		options.addOption("C", "cut", true,  "mutation count cutoff ");
		HelpFormatter formatter = new HelpFormatter();
		CommandLineParser parser = new BasicParser();
		CommandLine commandLine;
		try{
			commandLine = parser.parse(options, args);
		}catch (Exception e) {
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options] geneCountFile geneSetFile", options);
			return ;
		}
		List <String> argList = commandLine.getArgList();
		if(argList.size() != 2){
			formatter.printHelp(Reflection.getCallerClass( 1 ).getName() + " [options]  geneCountFile geneSetFile", options);
			return;
		}
		HotGenesetTest H = new HotGenesetTest();
		H.setGeneSet(argList.get(1));
		H.setMutationCount(argList.get(0));
		if(commandLine.hasOption("p")){
			H.setPcutoff(Integer.valueOf(commandLine.getOptionValue("p")));
		}
		
		if(commandLine.hasOption("c")){
			H.applyCeiling2MutationCount(Integer.valueOf(commandLine.getOptionValue("c")));
		}
		if(commandLine.hasOption("C")){
			H.cutSmallMutationCount(Integer.valueOf(commandLine.getOptionValue("C")));
		}
		H.perform();
		if(commandLine.hasOption("s")){
			H.printResultSimple();
		}else{
			H.printResult();
		}
	}
	
}
	